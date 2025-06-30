#!/bin/bash

# 查找 CPU 占用最大的进程
pid=$(ps -eo pid,%cpu,cmd --sort=-%cpu | head -n 2 | tail -n 1 | awk '{print $1}')

# 关闭进程
kill -9 $pid
echo "Process $pid is closed."
