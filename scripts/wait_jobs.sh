#!/bin/bash
# usage: wait_jobs.sh JOBID [JOBID...]  -- returns when none of the jobs is queued/running on the HPCC
ids=$(IFS=,; echo "$*")
while true; do
  n=$(ssh -o BatchMode=yes -o ConnectTimeout=20 hpcc2 "bash -lc 'squeue -h -u jstajich -j $ids 2>/dev/null | wc -l'" 2>/dev/null | tail -1)
  [ "$n" = "0" ] && { echo "jobs $ids finished $(date)"; exit 0; }
  sleep 120
done
