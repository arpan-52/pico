#!/bin/bash
#
# PICO v1 - Prompt Imaging of Cosmic Outbursts
# Trigger: parse burst config, drop a .job file, ensure the monitor is running.
# The monitor picks the job up and fires pico_image on the freest node.
#
# Usage:
#   pico.sh -c svfits_par.txt /lustre_scratch/.../RawVisi/DT_...
#
# Imaging params are NOT set here -- the node-side run job derives them from the
# observing band and appends them to the par file before calling pico_image.
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
MONITOR_SCRIPT="$SCRIPT_DIR/pico_monitor.py"
JOBS_PENDING_DIR="$SCRIPT_DIR/jobs_pending"

# Environment used to launch the monitor daemon (cluster layout, hardcoded).
TEST_IMAGING="/lustre_archive/apps/astrosoft/test_imaging"
START_ENV="$TEST_IMAGING/start.sh"
CONDA_ENV="pico"

log_message() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

show_usage() {
    cat << EOF
Usage: $0 -c config_file output_directory

Example:
    $0 -c svfits_par.txt /lustre_scratch/spotlight/data/SPLT_.../RawVisi/DT_...

Pipeline: pico_image (raw visibilities -> FITS) -> PNG
EOF
}

activate_env() {
    # Bring up the 'pico' micromamba env before anything runs (incl. our python3
    # calls). Hardcoded for the cluster; shell hook makes activate work non-login.
    source "$START_ENV"
    eval "$(micromamba shell hook --shell bash 2>/dev/null)" || true
    micromamba activate "$CONDA_ENV"
}

check_monitor() { pgrep -f "$MONITOR_SCRIPT" > /dev/null 2>&1; }

start_monitor() {
    log_message "Starting PICO monitor daemon..."
    [ -f "$MONITOR_SCRIPT" ] || { log_message "ERROR: monitor not found: $MONITOR_SCRIPT"; exit 1; }
    nohup bash -c "
source '$START_ENV'
eval \"\$(micromamba shell hook --shell bash 2>/dev/null)\" || true
micromamba activate $CONDA_ENV
python3 '$MONITOR_SCRIPT'
" > "$SCRIPT_DIR/monitor_startup.log" 2>&1 &
    sleep 3
    if check_monitor; then
        log_message "Monitor daemon started."
    else
        log_message "ERROR: monitor failed to start (see monitor_startup.log)"
        exit 1
    fi
}

extract_params() {
    local config_file="$1"
    [ -f "$config_file" ] || { log_message "ERROR: config not found: $config_file"; exit 1; }

    BURST_MJD=$(grep -Po 'BURST_MJD\s+\K[\d.]+' "$config_file" | head -1)
    BURST_RA=$(grep -Po 'BURST_RA\s+\K[\d.-]+' "$config_file" | head -1)
    BURST_DEC=$(grep -Po 'BURST_DEC\s+\K[\d.-]+' "$config_file" | head -1)

    # Observing band centre (Hz) from FREQ_SET = freq0:freq1:nchan
    local freq_set f0 f1
    freq_set=$(grep -Po 'FREQ_SET\s+\K[\d.:eE+-]+' "$config_file" | head -1)
    if [ -n "$freq_set" ]; then
        f0=${freq_set%%:*}
        f1=$(echo "$freq_set" | cut -d: -f2)
        BURST_FREQ=$(python3 -c "print(int((float('$f0')+float('$f1'))/2/1e6))")
    else
        BURST_FREQ=1400
        log_message "WARNING: FREQ_SET not found, defaulting band centre to 1400 MHz"
    fi

    [ -n "$BURST_MJD" ] || { log_message "ERROR: BURST_MJD not found in config"; exit 1; }
    # burst_<int>_<frac>
    BURST_NAME=$(python3 -c "p='$BURST_MJD'.split('.'); print('burst_'+p[0]+'_'+(p[1] if len(p)>1 else '0'))")

    log_message "Parsed: MJD=$BURST_MJD NAME=$BURST_NAME RA=$BURST_RA DEC=$BURST_DEC band_centre=${BURST_FREQ}MHz"
}

submit_job() {
    local config_file="$1" output_dir="$2"
    local job_id="${BURST_NAME}_$(date +%s)"

    mkdir -p "$JOBS_PENDING_DIR" "$output_dir"
    local job_log="$output_dir/pico_${job_id}.log"
    {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUBMIT] Job ID: $job_id"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUBMIT] Config: $(realpath "$config_file")"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUBMIT] Output: $(realpath "$output_dir")"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUBMIT] Burst MJD: $BURST_MJD  band_centre=${BURST_FREQ}MHz"
    } > "$job_log"

    cat > "$JOBS_PENDING_DIR/${job_id}.job" << EOF
JOB_ID=$job_id
BURST_NAME=$BURST_NAME
CONFIG_FILE=$(realpath "$config_file")
OUTPUT_DIR=$(realpath "$output_dir")
JOB_LOG=$job_log
BURST_MJD=$BURST_MJD
BURST_RA=$BURST_RA
BURST_DEC=$BURST_DEC
BURST_FREQ=$BURST_FREQ
SUBMIT_TIME=$(date '+%Y-%m-%d %H:%M:%S')
SUBMIT_HOST=$(hostname)
SUBMIT_USER=$(whoami)
EOF

    log_message "Job submitted: $job_id"
    log_message "  Job log: $job_log"
    echo ""
    echo "Watch:  tail -f $job_log"
    echo "Status: python3 $SCRIPT_DIR/pico_status.py pressure"
    echo ""
}

main() {
    local config_file="" output_dir=""
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--config) config_file="$2"; shift 2 ;;
            -h|--help)   show_usage; exit 0 ;;
            *)
                if [ -z "$output_dir" ]; then output_dir="$1"; shift
                else echo "ERROR: unknown arg: $1"; show_usage; exit 1; fi ;;
        esac
    done
    [ -n "$config_file" ] && [ -n "$output_dir" ] || { echo "ERROR: need -c config and output dir"; show_usage; exit 1; }

    mkdir -p "$JOBS_PENDING_DIR"
    activate_env
    extract_params "$config_file"
    check_monitor || start_monitor
    submit_job "$config_file" "$output_dir"
}

main "$@"
