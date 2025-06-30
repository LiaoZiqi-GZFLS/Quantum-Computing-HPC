#!/bin/bash
#SBATCH --job-name=sustechhpc-simulate
#SBATCH --partition=8175m          # 提交到 CPU 队列（partition 名）
#SBATCH --nodes=1                  # 所需节点数
#SBATCH --ntasks=1                 # 总任务数（通常 = 核心数）
#SBATCH --cpus-per-task=48         # 每个任务使用的 CPU 核心数
#SBATCH --time=10:00:00            # 最长运行时间（格式：hh:mm:ss）
#SBATCH --output=/work/sustcsc_02/workspace/output/slurm_%j.out      # 标准输出/错误日志（%j 会替换为作业ID）
#SBATCH --export=NONE              # 不继承环境变量

source /work/share/intel/oneapi-2023.1.0/setvars.sh
/work/sustcsc_02/workspace/simulate /work/share/simulate/input.bin

