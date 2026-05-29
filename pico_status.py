#!/usr/bin/env python3
"""
PICO v1 status / pressure viewer.

  pico_status.py pressure          # per-node load (active jobs) + queue summary
  pico_status.py pressure --watch  # refresh every 2s
  pico_status.py overview          # job status counts
  pico_status.py jobs [--status running] [--limit 20]
  pico_status.py detail <job_id>
  pico_status.py logs <job_id> [--lines 50]
"""

import os
import sys
import time
import sqlite3
import argparse
from pathlib import Path

DB = Path(__file__).resolve().parent / "pico_jobs.db"


def conn():
    if not DB.exists():
        sys.exit(f"Database not found: {DB} (has the monitor run yet?)")
    c = sqlite3.connect(str(DB), timeout=30)
    c.row_factory = sqlite3.Row
    return c


def status_counts(c):
    rows = c.execute("SELECT status, COUNT(*) n FROM jobs GROUP BY status").fetchall()
    return {r["status"]: r["n"] for r in rows}


def cmd_pressure(args):
    def render():
        with conn() as c:
            srv = c.execute('''SELECT server, active_jobs, total_done, last_updated
                FROM server_status ORDER BY server ASC''').fetchall()
            counts = status_counts(c)
        out = []
        out.append("=== Node pressure ===")
        out.append(f"{'Node':<10}{'Active':>7}{'Done':>7}   {'Last update':<19}  Load")
        out.append("-" * 58)
        peak = max((r["active_jobs"] for r in srv), default=0) or 1
        for r in srv:
            bar = "#" * int(round(10 * r["active_jobs"] / peak)) if r["active_jobs"] else ""
            lu = (r["last_updated"] or "")[:19]
            out.append(f"{r['server']:<10}{r['active_jobs']:>7}{r['total_done']:>7}   {lu:<19}  {bar}")
        out.append("-" * 58)
        total_active = sum(r["active_jobs"] for r in srv)
        out.append(f"Total active: {total_active}    "
                   f"PENDING={counts.get('PENDING',0)} RUNNING={counts.get('RUNNING',0)} "
                   f"COMPLETED={counts.get('COMPLETED',0)} FAILED={counts.get('FAILED',0)}")
        return "\n".join(out)

    if getattr(args, "watch", False):
        try:
            while True:
                os.system("clear")
                print(render())
                print("\n(ctrl-c to exit)")
                time.sleep(2)
        except KeyboardInterrupt:
            pass
    else:
        print(render())


def cmd_overview(args):
    with conn() as c:
        counts = status_counts(c)
    print("=== PICO v1 overview ===")
    for s in ["PENDING", "RUNNING", "COMPLETED", "FAILED"]:
        print(f"  {s:<10}: {counts.get(s, 0):4d}")


def cmd_jobs(args):
    q = "SELECT job_id, server, status, submit_time, burst_mjd FROM jobs"
    p = []
    if args.status:
        q += " WHERE status=?"
        p.append(args.status.upper())
    q += " ORDER BY submit_time DESC"
    if args.limit:
        q += f" LIMIT {int(args.limit)}"
    with conn() as c:
        rows = c.execute(q, p).fetchall()
    if not rows:
        print("No jobs found")
        return
    print(f"{'Job ID':<34}{'Node':<10}{'Status':<11}{'Submitted':<20}")
    print("-" * 75)
    for r in rows:
        print(f"{r['job_id'][:33]:<34}{(r['server'] or 'N/A'):<10}"
              f"{r['status']:<11}{(r['submit_time'] or 'N/A'):<20}")


def cmd_detail(args):
    with conn() as c:
        r = c.execute("SELECT * FROM jobs WHERE job_id=?", (args.job_id,)).fetchone()
    if not r:
        print(f"Job {args.job_id} not found")
        return
    for k in r.keys():
        print(f"{k:<14}: {r[k]}")
    log = r["job_log"]
    if log and Path(log).exists():
        print("\n=== last 20 log lines ===")
        print("".join(open(log).readlines()[-20:]).rstrip())


def cmd_logs(args):
    with conn() as c:
        r = c.execute("SELECT job_log FROM jobs WHERE job_id=?", (args.job_id,)).fetchone()
    if not r or not r["job_log"] or not Path(r["job_log"]).exists():
        print("Log not found")
        return
    print("".join(open(r["job_log"]).readlines()[-args.lines:]).rstrip())


def main():
    ap = argparse.ArgumentParser(description="PICO v1 status / pressure viewer")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("pressure"); p.add_argument("--watch", action="store_true")
    sub.add_parser("overview")
    j = sub.add_parser("jobs")
    j.add_argument("--status", choices=["pending", "running", "completed", "failed"])
    j.add_argument("--limit", type=int, default=20)
    d = sub.add_parser("detail"); d.add_argument("job_id")
    l = sub.add_parser("logs"); l.add_argument("job_id"); l.add_argument("--lines", type=int, default=50)

    args = ap.parse_args()
    {
        None: cmd_pressure, "pressure": cmd_pressure, "overview": cmd_overview,
        "jobs": cmd_jobs, "detail": cmd_detail, "logs": cmd_logs,
    }[args.cmd](args)


if __name__ == "__main__":
    main()
