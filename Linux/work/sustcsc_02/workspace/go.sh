source /work/share/intel/oneapi-2023.1.0/setvars.sh
icpx -std=c++17 -xHost -qopenmp -O3 /work/sustcsc_02/workspace/simulate.cpp /work/sustcsc_02/workspace/driver.o -o /work/sustcsc_02/workspace/simulate
