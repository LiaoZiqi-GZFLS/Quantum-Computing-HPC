#!/bin/bash

# 提交作业并获取作业 ID
JOB_ID=$(sbatch /work/sustcsc_02/workspace/job.sh | awk '{print $4}')
echo "Submitted job with ID: $JOB_ID"

# 等待作业结束
while squeue -j $JOB_ID | grep -q "$JOB_ID"; do
  squeue -u sustcsc_02
  sleep 10
done

# 输出 .out 文件内容
OUT_FILE="/work/sustcsc_02/workspace/output/slurm_${JOB_ID}.out"
if [ -f "$OUT_FILE" ]; then
  echo "Content of $OUT_FILE:"
  cat "$OUT_FILE"
else
  echo "Output file $OUT_FILE not found."
fi

sacct -j $JOB_ID
