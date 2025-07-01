#include <chrono>
#include <dlfcn.h>
#include <cstdio> // 包含 cstdio 头文件
#include <stdlib.h>

using namespace std::chrono;

// 原始的 std::chrono::steady_clock::now 函数指针
using original_now_t = steady_clock::time_point(*)();
original_now_t original_now = nullptr;

// 拦截后的 std::chrono::steady_clock::now 函数
steady_clock::time_point now() {
    if (getenv("LD_PRELOAD") == NULL) {
        return steady_clock::time_point();
    }
    unsetenv("LD_PRELOAD");
    if (!original_now) {
        // 动态加载 libstdc++.so.6 并获取原始的 std::chrono::steady_clock::now 函数
        void* handle = dlopen("libstdc++.so.6", RTLD_LAZY);
        if (!handle) {
            // 使用 fprintf 而不是 std::fprintf
            fprintf(stderr, "Failed to load libstdc++.so.6: %s\n", dlerror());
            return steady_clock::time_point();
        }

        original_now = (original_now_t)dlsym(handle, "std::chrono::steady_clock::now");
        if (!original_now) {
            // 使用 fprintf 而不是 std::fprintf
            fprintf(stderr, "Failed to find std::chrono::steady_clock::now in libstdc++.so.6: %s\n", dlerror());
            return steady_clock::time_point();
        }
    }

    // 调用原始的 std::chrono::steady_clock::now 函数
    auto original_time = original_now();

    // 添加自定义延迟（例如，增加 100 毫秒）
    return original_time + std::chrono::milliseconds(100);
}