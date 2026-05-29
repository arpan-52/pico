#!/usr/bin/env python3
"""
PICO v1 Monitor - distributes pico_image jobs across the compute nodes.

Scans jobs_pending/ for .job files, records them in pico_jobs.db, then dispatches
each to the *freest* node (fewest active jobs) over ssh. Node load is tracked in
the server_status table: a job atomically increments its node's counter at
dispatch and decrements it on completion, so selection actually spreads load
across all nodes instead of always picking the first one.
"""

import os
import sys
import time
import sqlite3
import subprocess
import signal
import logging
import threading
from datetime import datetime
from pathlib import Path

SERVERS = ["rggpur00", "rggpur01", "rggpur02", "rggpur03",
           "rggpur04", "rggpur05", "rggpur06", "rggpur07"]


class PICOMonitor:
    def __init__(self, script_dir):
        self.script_dir = Path(script_dir)
        self.db_path = self.script_dir / "pico_jobs.db"
        self.jobs_pending_dir = self.script_dir / "jobs_pending"
        self.run_job_script = self.script_dir / "pico_run_job.sh"
        self.pico_image_bin = (self.script_dir / "pico_backend" /
                               "pico_image" / "build" / "pico_image")
        self.jobs_pending_dir.mkdir(exist_ok=True)

        self.servers = SERVERS
        self.active_jobs = {}
        self.shutdown = False
        # Serialises node selection + counter increment so two jobs dispatched in
        # the same tick don't both land on the same "freest" node.
        self.sched_lock = threading.Lock()

        self.setup_logging()
        self.init_database()

        signal.signal(signal.SIGTERM, self.signal_handler)
        signal.signal(signal.SIGINT, self.signal_handler)
        self.logger.info("PICO v1 monitor started in %s", self.script_dir)
        self.logger.info("pico_image binary: %s", self.pico_image_bin)

    def setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s - %(levelname)s - %(message)s",
            handlers=[logging.FileHandler(self.script_dir / "monitor.log"),
                      logging.StreamHandler()],
        )
        self.logger = logging.getLogger("pico_monitor")

    def signal_handler(self, signum, frame):
        self.logger.info("Received signal %s, shutting down...", signum)
        self.shutdown = True

    def db(self):
        conn = sqlite3.connect(str(self.db_path), timeout=30)
        conn.execute("PRAGMA busy_timeout=30000")
        return conn

    def init_database(self):
        with self.db() as conn:
            conn.execute('''
                CREATE TABLE IF NOT EXISTS jobs (
                    job_id TEXT PRIMARY KEY,
                    server TEXT,
                    output_dir TEXT,
                    job_log TEXT,
                    config_file TEXT,
                    status TEXT DEFAULT 'PENDING',
                    submit_time TEXT,
                    start_time TEXT,
                    end_time TEXT,
                    burst_mjd REAL,
                    burst_freq INTEGER,
                    error_msg TEXT
                )''')
            conn.execute('''
                CREATE TABLE IF NOT EXISTS server_status (
                    server TEXT PRIMARY KEY,
                    active_jobs INTEGER DEFAULT 0,
                    total_done INTEGER DEFAULT 0,
                    last_updated TEXT
                )''')
            for srv in self.servers:
                conn.execute('''INSERT OR IGNORE INTO server_status
                    (server, active_jobs, total_done, last_updated)
                    VALUES (?, 0, 0, ?)''', (srv, datetime.now().isoformat()))
            # Reset stale counters from a previous (crashed) run.
            conn.execute("UPDATE server_status SET active_jobs=0")
            conn.commit()
        self.logger.info("Database initialised")

    def scan_pending_jobs(self):
        for job_file in self.jobs_pending_dir.glob("*.job"):
            try:
                params = {}
                with open(job_file) as f:
                    for line in f:
                        line = line.strip()
                        if "=" in line:
                            k, v = line.split("=", 1)
                            params[k] = v
                job_id = params.get("JOB_ID")
                if not job_id:
                    continue
                with self.db() as conn:
                    exists = conn.execute(
                        "SELECT 1 FROM jobs WHERE job_id=?", (job_id,)).fetchone()
                    if not exists:
                        conn.execute('''INSERT INTO jobs
                            (job_id, output_dir, job_log, config_file, status,
                             submit_time, burst_mjd, burst_freq)
                            VALUES (?,?,?,?,'PENDING',?,?,?)''',
                            (job_id, params.get("OUTPUT_DIR"), params.get("JOB_LOG"),
                             params.get("CONFIG_FILE"), params.get("SUBMIT_TIME"),
                             float(params.get("BURST_MJD", 0)),
                             int(params.get("BURST_FREQ", 0))))
                        conn.commit()
                job_file.unlink()
                self.logger.info("Queued job: %s", job_id)
            except Exception as e:
                self.logger.error("Error reading %s: %s", job_file.name, e)

    def claim_freest_server(self, job_id):
        """Pick the node with the fewest active jobs and increment it atomically."""
        with self.sched_lock, self.db() as conn:
            row = conn.execute('''SELECT server FROM server_status
                ORDER BY active_jobs ASC, server ASC LIMIT 1''').fetchone()
            server = row[0] if row else self.servers[0]
            conn.execute('''UPDATE server_status
                SET active_jobs = active_jobs + 1, last_updated=?
                WHERE server=?''', (datetime.now().isoformat(), server))
            conn.commit()
        return server

    def release_server(self, server, ok):
        with self.db() as conn:
            conn.execute('''UPDATE server_status
                SET active_jobs = MAX(0, active_jobs - 1),
                    total_done = total_done + ?, last_updated=?
                WHERE server=?''',
                (1 if ok else 0, datetime.now().isoformat(), server))
            conn.commit()

    def get_pending_jobs(self):
        with self.db() as conn:
            return conn.execute('''SELECT job_id, output_dir, job_log, config_file, burst_mjd
                FROM jobs WHERE status='PENDING' ORDER BY submit_time ASC''').fetchall()

    def execute_job(self, job_data):
        job_id, output_dir, job_log, config_file, burst_mjd = job_data
        server = self.claim_freest_server(job_id)
        self.logger.info("Dispatching %s -> %s", job_id, server)

        with self.db() as conn:
            conn.execute('''UPDATE jobs SET server=?, status='RUNNING', start_time=?
                WHERE job_id=?''', (server, datetime.now().isoformat(), job_id))
            conn.commit()
        with open(job_log, "a") as f:
            f.write(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] [MONITOR] Assigned to {server}\n")

        ok = False
        try:
            cmd = " ".join([
                "ssh", server, f"bash {self.run_job_script}",
                job_id, output_dir, config_file, str(burst_mjd),
                job_log, str(self.pico_image_bin),
            ])
            self.logger.info("Running: %s", cmd)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True)
            stdout, _ = proc.communicate()
            if stdout:
                with open(job_log, "a") as f:
                    f.write(stdout)
            ok = proc.returncode == 0
            status = "COMPLETED" if ok else "FAILED"
            err = None if ok else f"run_job exit {proc.returncode}"
            with self.db() as conn:
                conn.execute('''UPDATE jobs SET status=?, end_time=?, error_msg=?
                    WHERE job_id=?''', (status, datetime.now().isoformat(), err, job_id))
                conn.commit()
            (self.logger.info if ok else self.logger.error)(
                "Job %s %s on %s", job_id, status, server)
        except Exception as e:
            self.logger.error("Exception running %s: %s", job_id, e)
            with self.db() as conn:
                conn.execute('''UPDATE jobs SET status='FAILED', end_time=?, error_msg=?
                    WHERE job_id=?''', (datetime.now().isoformat(), str(e), job_id))
                conn.commit()
        finally:
            self.release_server(server, ok)
            self.active_jobs.pop(job_id, None)

    def run(self):
        self.logger.info("Monitor PID: %s", os.getpid())
        while not self.shutdown:
            try:
                self.scan_pending_jobs()
                for job_data in self.get_pending_jobs():
                    if self.shutdown:
                        break
                    job_id = job_data[0]
                    if job_id in self.active_jobs:
                        continue
                    t = threading.Thread(target=self.execute_job,
                                         args=(job_data,), daemon=True)
                    self.active_jobs[job_id] = t
                    t.start()
                    self.logger.info("Launched: %s", job_id)
                time.sleep(1)
            except Exception as e:
                self.logger.error("Main loop error: %s", e)
                time.sleep(5)
        self.logger.info("Monitor shutting down")


def main():
    monitor = PICOMonitor(os.path.dirname(os.path.abspath(__file__)))
    try:
        monitor.run()
    except KeyboardInterrupt:
        monitor.logger.info("Stopped by user")
    except Exception as e:
        monitor.logger.error("Monitor crashed: %s", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
