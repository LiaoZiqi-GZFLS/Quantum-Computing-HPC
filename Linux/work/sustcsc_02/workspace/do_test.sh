#!/bin/sh

# 初始化计数器
correct=0
total=20

for i in $(seq 1 $total)
do
	echo "times: $i"

	# 运行两个脚本并将输出写入文件
	/work/sustcsc_02/workspace/gen $i /work/sustcsc_02/workspace/input.bin
	/work/sustcsc_02/workspace/run.sh > /work/sustcsc_02/workspace/result.txt
	/work/sustcsc_02/workspace/ref_run.sh > /work/sustcsc_02/workspace/ref_result.txt

	# 获取两个文件的第一行
	line1=$(head -n 1 /work/sustcsc_02/workspace/result.txt)
	line2=$(head -n 1 /work/sustcsc_02/workspace/ref_result.txt)
	
	# 比较两个文件的第一行
	if [ "$line1" = "$line2" ]; then
        	((correct++))
	else
		echo $line1$
		echo $line2$
	fi
done

# 计算准确率
accuracy=$(echo "scale=2; $correct / $total * 100" | bc)
echo "Accuracy: $accuracy%"
