#!/bin/bash

# 查找 CPU 占用最大的进程
pid=$(ps -eo pid,%cpu,cmd --sort=-%cpu | head -n 2 | tail -n 1 | awk '{print $1}')

# 输出进程信息
echo "CPU 占用最大的进程 PID: $pid"

# 确认是否要关闭该进程
read -p "是否要关闭该进程? (y/n): " confirm

if [[ "$confirm" == "y" || "$confirm" == "Y" ]]; then
    # 关闭进程
    kill -9 $pid
    echo "进程 $pid 已被关闭。"
else
    echo "取消操作。"
fi
