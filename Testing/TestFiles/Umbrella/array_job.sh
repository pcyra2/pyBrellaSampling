#!/bin/bash
RUNLINE=$(cat $ARRAY_TASKFILE | head -n $(($SLURM_ARRAY_TASK_ID*1)) | tail -n 1)
eval "$RUNLINE wait"


