#include <iostream>
#include <cstdlib>
#include <Windows.h>
#include <iomanip>

// 方法1：列优先遍历（缓存不友好）
void method1(double** matrix, double* vector, double* result, int n) {
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += matrix[i][j] * vector[i]; // 列优先访问
        }
        result[j] = sum;
    }
}

// 方法2：显式列优先累加
void method2(double** b, double* a, double* sum, int n) {
    for (int i = 0; i < n; ++i) sum[i] = 0.0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            sum[i] += b[j][i] * a[j]; // 列优先访问
        }
    }
}

// 方法3：循环展开
void method3(double** b, double* a, double* sum, int n) {
    for (int i = 0; i < n; ++i) sum[i] = 0.0;
    for (int j = 0; j < n; ++j) {
    int i;
    // 主循环，每次处理4个元素
    for (i = 0; i <= n - 4; i += 4) {
         sum[i]     += b[j][i]     * a[j];
         sum[i + 1] += b[j][i + 1] * a[j];
         sum[i + 2] += b[j][i + 2] * a[j];
         sum[i + 3] += b[j][i + 3] * a[j];
    }
    // 处理剩余不足4个的元素
    for (; i < n; ++i) {
         sum[i] += b[j][i] * a[j];
    }
}

}

int main() {
    const int n = 20000;      // 矩阵维度（测试大尺寸最坏情况）
    const int trials = 1;   // 测试次数

    // 动态分配内存（模拟列优先存储）
    double** matrix = new double*[n]; // 方法1的矩阵（行优先存储）
    double** b = new double*[n];      // 方法2的矩阵（列优先存储）
    double* vector = new double[n];   // 输入向量
    double* a = new double[n];        // 输入向量副本
    double* result1 = new double[n];  // 方法1结果
    double* result2 = new double[n];  // 方法2结果
    double* result3 = new double[n];  // 方法3结果

    // 分别为 matrix 和 b 分配每一行的内存
    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[n];
    }
    for (int i = 0; i < n; ++i) {
        b[i] = new double[n];
    }

    // 初始化数据
    for (int i = 0; i < n; ++i) {
        vector[i] = static_cast<double>(rand()) / RAND_MAX;
        a[i] = vector[i];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = static_cast<double>(rand()) / RAND_MAX;
            b[j][i] = matrix[i][j]; // 将 b 设为 matrix 的转置（模拟列优先）
        }
    }

    // 高精度计时函数
    long long head1, tail1, head2, tail2, freq, head3, tail3;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    // 测试方法1
    double time1 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head1);
        method1(matrix, vector, result1, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);
        time1 += (tail1 - head1) *1000 / freq;
    }
    // 测试方法2
    double time2 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head2);
        method2(b, a, result2, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
        time2 += (tail2 - head2) *1000 / freq;
    }

    // 测试方法3
    double time3 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head3);
        method3(b, a, result3, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail3);
        time3 += (tail3 - head3) *1000 / freq;
    }

    // 输出结果
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Method1 (列优先访问) 平均时间: " << time1 / trials << " 毫秒" << std::endl;
    std::cout << "Method2 (显式列累加) 平均时间: " << time2 / trials << " 毫秒" << std::endl;
    std::cout << "Method3 (循环优化) 平均时间: " << time3 / trials << " 毫秒" << std::endl;

    // 释放内存
    for (int i = 0; i < n; ++i) {
        delete[] matrix[i];
        delete[] b[i];
    }
    delete[] matrix;
    delete[] b;
    delete[] vector;
    delete[] a;
    delete[] result1;
    delete[] result2;
    delete[] result3;

    return 0;
}
