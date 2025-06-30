#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <complex>
#include <memory>
#include <cstring>

using namespace std;

extern "C" {
    void __printf_chk(double, int, const char*, ...);
    void __fprintf_chk(FILE*, int, const char*, ...);
    void std::__throw_length_error(const char* msg);
    void operator_delete(void* ptr);
    void* operator_new(size_t size);
}

int main(int argc, char** argv) {
    char* alpha_real_ptr;
    char* alpha_imag_ptr;
    char* beta_real_ptr;
    char* beta_imag_ptr;
    FILE* file_stream;
    char* data_ptr;
    size_t data_size;
    long start_time;
    long end_time;
    long fs_offset;

    // 初始化栈保护变量
    fs_offset = *(long*)(in_FS_OFFSET + 0x28);

    // 检查参数个数
    if (argc == 2) {
        // 打开文件
        file_stream = fopen(argv[1], "rb");
        if (file_stream == nullptr) {
            perror("Failed to open file");
            goto error_exit;
        }

        // 读取数据大小
        fread(&data_size, sizeof(data_size), 1, file_stream);

        // 检查数据大小
        if (data_size > (size_t)-1 / sizeof(double)) {
            if (fs_offset == *(long*)(in_FS_OFFSET + 0x28)) {
                std::__throw_length_error("cannot create std::vector larger than max_size()");
            }
            goto stack_check_fail;
        }

        // 分配内存
        if (data_size == 0) {
            data_ptr = nullptr;
            alpha_real_ptr = nullptr;
            alpha_imag_ptr = nullptr;
            beta_real_ptr = nullptr;
            beta_imag_ptr = nullptr;
        } else {
            data_ptr = static_cast<char*>(operator_new(data_size));
            memset(data_ptr, 0, data_size);
            alpha_real_ptr = data_ptr;
            alpha_imag_ptr = data_ptr + sizeof(double);
            beta_real_ptr = data_ptr + 2 * sizeof(double);
            beta_imag_ptr = data_ptr + 3 * sizeof(double);
        }

        // 尝试读取数据
        try {
            fread(data_ptr, 1, data_size, file_stream);
        } catch (...) {
            // 捕获异常并处理
            fclose(file_stream);
            operator_delete(data_ptr);
            goto error_exit;
        }

        fclose(file_stream);

        // 获取当前时间
        start_time = std::chrono::_V2::system_clock::now();

        // 调用模拟函数
        simulate(data_size, data_ptr, reinterpret_cast<complex*>(alpha_real_ptr), reinterpret_cast<complex*>(beta_real_ptr));

        end_time = std::chrono::_V2::system_clock::now();

        // 打印结果
        __printf_chk(2, "Final state: alpha = %.12f + %.12fi, beta = %.12f + %.12fi\n",
                     *reinterpret_cast<double*>(alpha_real_ptr),
                     *reinterpret_cast<double*>(alpha_imag_ptr),
                     *reinterpret_cast<double*>(beta_real_ptr),
                     *reinterpret_cast<double*>(beta_imag_ptr));

        // 打印时间
        __printf_chk(2, "Time taken: %.2f ms\n", ((double)(end_time - start_time) / 1000000.0));

        // 释放内存
        operator_delete(data_ptr);

        return 0;
    } else {
        // 打印用法信息
        __fprintf_chk(stderr, 2, "Usage: %s <input_file>\n", argv[0]);
        goto error_exit;
    }

error_exit:
    return 1;

stack_check_fail:
    __stack_chk_fail();
    return 1;
}