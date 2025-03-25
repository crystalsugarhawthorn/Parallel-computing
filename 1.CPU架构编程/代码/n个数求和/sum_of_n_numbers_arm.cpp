#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <iomanip>
#include <sys/time.h>  // 替换 Windows 计时函数
#include <locale.h>    // 支持中文字符

using namespace std;

// 简单求和
void simple_sum(double a[], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += a[i];
    }
}

// 递归求和
void recursive_sum(double a[], int n) {
    if (n == 1) return;
    for (int i = 0; i < n/2; ++i) {
        a[i] += a[n - i - 1];
    }
    recursive_sum(a, n/2);
}

// 迭代求和
void iterative_sum(double a[], int n) {
    for (int m = n; m > 1; m /= 2) {
        for (int i = 0; i < m/2; ++i) {
            a[i] = a[2*i] + a[2*i + 1];
        }
    }
}

// 循环展开求和
void iterativePlus_sum(double a[], int n) {
    double sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
    int i;
    for (i = 0; i <= n - 4; i += 4) {
        sum0 += a[i];
        sum1 += a[i + 1];
        sum2 += a[i + 2];
        sum3 += a[i + 3];
    }
    for (; i < n; ++i) {
        sum0 += a[i];
    }
    double sum = sum0 + sum1 + sum2 + sum3;
}

// 高精度计时函数 (Linux/POSIX)
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000.0 + tv.tv_usec / 1000.0;
}

int main() {
    setlocale(LC_ALL, ""); // 支持中文字符输出

    const int n = 10000000;   // 数组大小
    const int trials = 10;      // 测试次数

    // 初始化数据（确保每次测试数据一致）
    vector<double> arr(n);
    srand(time(NULL));
    for (int i = 0; i < n; ++i) {
        arr[i] = static_cast<double>(rand()) / RAND_MAX;
    }

    // 为每个方法分配独立的数据副本
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *d = new double[n];
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = c[i] = d[i] = arr[i];
    }

    // 测试方法1：简单求和
    double start1 = get_time();
    for (int t = 0; t < trials; ++t) {
        simple_sum(a, n);
    }
    double time1 = get_time() - start1;

    // 测试方法2：递归求和
    double start2 = get_time();
    for (int t = 0; t < trials; ++t) {
        recursive_sum(b, n);
    }
    double time2 = get_time() - start2;

    // 测试方法3：迭代求和
    double start3 = get_time();
    for (int t = 0; t < trials; ++t) {
        iterative_sum(c, n);
    }
    double time3 = get_time() - start3;

    // 测试方法4：循环展开求和
    double start4 = get_time();
    for (int t = 0; t < trials; ++t) {
        iterativePlus_sum(d, n);
    }
    double time4 = get_time() - start4;

    // 输出结果
    cout << fixed << setprecision(6);
    cout << "性能测试 (n=" << n << ", trials=" << trials << "):\n";
    cout << "1. 简单求和耗时: " << time1 / trials << " ms\n";
    cout << "2. 递归求和耗时: " << time2 / trials << " ms\n";
    cout << "3. 迭代求和耗时: " << time3 / trials << " ms\n";
    cout << "4. 循环展开求和耗时: " << time4 / trials << " ms\n";

    // 释放内存
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    return 0;
}