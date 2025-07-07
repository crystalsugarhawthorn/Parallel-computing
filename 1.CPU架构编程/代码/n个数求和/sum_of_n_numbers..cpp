#include <iostream>
#include <cstdlib>
#include <vector>
#include <windows.h>
#include <iomanip>

using namespace std;

// 平凡算法
void simple_sum(double a[], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += a[i];
    }
}

// 递归算法
void recursive_sum(double a[], int n) {
    if (n == 1) return;           // 递归终止条件
    // 折叠数组：前n/2元素与后n/2元素相加
    for (int i = 0; i < n/2; ++i) {
        a[i] += a[n - i - 1];    // 镜像元素相加
    }
    recursive_sum(a, n/2);       // 递归处理前半段
}

// 二重循环
void iterative_sum(double a[], int n) {
    for (int m = n; m > 1; m /= 2) {       // 外层控制数组规模
        for (int i = 0; i < m/2; ++i) {    // 内层处理相邻元素
            a[i] = a[2*i] + a[2*i + 1];    // 相邻元素压缩存储
        }
    }
}

//循环展开
void iterativePlus_sum(double a[], int n) {
    int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
    int i;
    // 主循环，每次处理4个元素
    for (i = 0; i <= n - 4; i += 4) {
        sum0 += a[i];
        sum1 += a[i + 1];
        sum2 += a[i + 2];
        sum3 += a[i + 3];
    }
    // 处理剩余不足4个的元素
    for (; i < n; ++i) {
        sum0 += a[i];
    }
    // 最终结果为各累加器之和
    int sum = sum0 + sum1 + sum2 + sum3;
}


int main()
{
    const int n = 100000000;   // 测试数据规模
    const int trials = 1;   // 测试轮次

    // 生成随机测试数据 (预分配连续内存)
    vector<double> arr(n);
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *d = new double[n];
    for (int i = 0; i < n; ++i){
        arr[i] = static_cast<double>(rand()) / RAND_MAX;
        a[i] = b[i] = c[i] = d[i] = arr[i];
    }

    // 高精度计时函数
    long long head1, tail1, head2, tail2, head3, tail3, head4, tail4, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    // 测试方法1
    double time1 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head1);
        simple_sum(a, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);
        time1 += (tail1 - head1) *1000 / freq;
    }
    // 测试方法2
    double time2 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head2);
        recursive_sum(b, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
        time2 += (tail2 - head2) *1000 / freq;
    }
    // 测试方法3
    double time3 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head3);
        iterative_sum(c, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail3);
        time3 += (tail3 - head3) *1000 / freq;
    }
    // 测试方法3
    double time4 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head4);
        iterativePlus_sum(a, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail4);
        time4 += (tail4 - head4) *1000 / freq;
    }
    // 输出结果
    cout << fixed << setprecision(6);
    cout << "测试结果 (n=" << n << ", trials=" << trials << "):\n";
    cout << "1. 平凡算法平均耗时: " << time1 / trials << " ms\n";
    cout << "2. 递归算法平均耗时: " << time2 / trials << " ms\n";
    cout << "3. 二重循环平均耗时: " << time3 / trials << " ms\n";
    cout << "4. 循环展开平均耗时: " << time4 / trials << " ms\n";

    // 测试方法4
    for (int t = 0; t < trials; ++t) {
        iterativePlus_sum(d, n);
    }


    // 释放内存
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    return 0;
}
