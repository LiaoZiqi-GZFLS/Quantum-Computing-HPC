#!/bin/bash

#SBATCH --job-name=sustechhpc-simulate
#SBATCH --partition=8175m          # 提交到 CPU 队列（partition 名）
#SBATCH --nodes=1                  # 所需节点数
#SBATCH --ntasks=1                 # 总任务数（通常 = 核心数）
#SBATCH --cpus-per-task=48         # 每个任务使用的 CPU 核心数
#SBATCH --time=10:00:00            # 最长运行时间（格式：hh:mm:ss）
#SBATCH --output=/work/sustcsc_02/workspace/output/slurm_%j.out      # 标准输出/错误日志（%j 会替换为作业ID）
#SBATCH --export=NONE              # 不继承环境变量


# 获取当前时间的毫秒部分
milliseconds=$(date +'%3N')
# 取前 4 位毫秒数
four_digit_num=$(printf "%04d" $milliseconds)
echo $four_digit_num

# 要检测的目录路径
directory_path="/work/sustcsc_02/workspace/code"

# 初始化一个空变量来存储 .cpp 文件
cpp_files=""

# 遍历目录中的所有文件和文件夹
for file in "$directory_path"/*; do
    # 检查是否为文件（排除目录）
    if [ -f "$file" ]; then
        # 提取文件扩展名
        file_extension="${file##*.}"
        # 检查文件扩展名是否为 cpp
        if [ "$file_extension" = "cpp" ]; then
            # 如果变量是空的，直接赋值
            if [ -z "$cpp_files" ]; then
                cpp_files="$file"
		break
            # 如果变量已有值，追加空格和新的文件路径
            else
                break
            fi
        fi
    fi
done

# 输出结果（如果需要）
echo "Found .cpp files:"
echo "$cpp_files"

#go.sh
source /work/share/intel/oneapi-2023.1.0/setvars.sh
icpx -std=c++17 -xHost -qopenmp -O3 $cpp_files /work/sustcsc_02/workspace/driver.o -o /work/sustcsc_02/workspace/code/simulate_$four_digit_num
rm $cpp_files

#do_test.sh
# 初始化计数器
correct=0
total=20

for i in $(seq 1 $total)
do
	echo "times: $i"

	# 运行两个脚本并将输出写入文件
	/work/sustcsc_02/workspace/gen $i /work/sustcsc_02/workspace/code/input_$four_digit_num.bin
	#/work/sustcsc_02/workspace/run.sh > /work/sustcsc_02/workspace/result.txt
	/work/sustcsc_02/workspace/code/simulate_$four_digit_num /work/sustcsc_02/workspace/code/input_$four_digit_num.bin > /work/sustcsc_02/workspace/code/result_$four_digit_num.txt
	#/work/sustcsc_02/workspace/ref_run.sh > /work/sustcsc_02/workspace/ref_result.txt
	/work/sustcsc_02/workspace/ref_simulate /work/sustcsc_02/workspace/code/input_$four_digit_num.bin > /work/sustcsc_02/workspace/code/ref_result_$four_digit_num.txt
	
	# 获取两个文件的第一行
	line1=$(head -n 1 /work/sustcsc_02/workspace/code/result_$four_digit_num.txt)
	line2=$(head -n 1 /work/sustcsc_02/workspace/code/ref_result_$four_digit_num.txt)

	# 比较两个文件的第一行
	if [ "$line1" = "$line2" ]; then
        	((correct++))
	else
		echo $line1$
		echo $line2$
	fi
done

rm /work/sustcsc_02/workspace/code/result_$four_digit_num.txt
rm /work/sustcsc_02/workspace/code/ref_result_$four_digit_num.txt
rm /work/sustcsc_02/workspace/code/input_$four_digit_num.bin
rm /work/sustcsc_02/workspace/code/simulate_$four_digit_num

# 计算准确率
accuracy=$(echo "scale=2; $correct / $total * 100" | bc)
echo "Accuracy: $accuracy%"

