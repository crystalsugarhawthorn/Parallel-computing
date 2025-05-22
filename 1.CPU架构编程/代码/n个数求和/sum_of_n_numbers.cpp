#include <iostream>
#include <cstdlib>
#include <vector>
#include <windows.h>
#include <iomanip>

using namespace std;

// ?????
void simple_sum(double a[], int n) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += a[i];
    }
}

// ?????
void recursive_sum(double a[], int n) {
    if (n == 1) return;           // ??????????
    // ??????ï“?n/2??????n/2??????
    for (int i = 0; i < n/2; ++i) {
        a[i] += a[n - i - 1];    // ??????????
    }
    recursive_sum(a, n/2);       // ??ø„??????
}

// ???????
void iterative_sum(double a[], int n) {
    for (int m = n; m > 1; m /= 2) {       // ????????????
        for (int i = 0; i < m/2; ++i) {    // ????????????
            a[i] = a[2*i] + a[2*i + 1];    // ???????????›¥
        }
    }
}

//??????
void iterativePlus_sum(double a[], int n) {
    int sum0 = 0, sum1 = 0, sum2 = 0, sum3 = 0;
    int i;
    // ?????????¦Ä???4?????
    for (i = 0; i <= n - 4; i += 4) {
        sum0 += a[i];
        sum1 += a[i + 1];
        sum2 += a[i + 2];
        sum3 += a[i + 3];
    }
    // ?????????4???????
    for (; i < n; ++i) {
        sum0 += a[i];
    }
    // ?????????????????
    int sum = sum0 + sum1 + sum2 + sum3;
}


int main()
{
    const int n = 100000000;   // ??????????
    const int trials = 1;   // ???????

    // ??????????????? (????????????)
    vector<double> arr(n);
    double *a = new double[n];
    double *b = new double[n];
    double *c = new double[n];
    double *d = new double[n];
    for (int i = 0; i < n; ++i){
        arr[i] = static_cast<double>(rand()) / RAND_MAX;
        a[i] = b[i] = c[i] = d[i] = arr[i];
    }

    // ???????????
    long long head1, tail1, head2, tail2, head3, tail3, head4, tail4, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    // ???????1
    double time1 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head1);
        simple_sum(a, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);
        time1 += (tail1 - head1) *1000 / freq;
    }
    // ???????2
    double time2 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head2);
        recursive_sum(b, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
        time2 += (tail2 - head2) *1000 / freq;
    }
    // ???????3
    double time3 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head3);
        iterative_sum(c, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail3);
        time3 += (tail3 - head3) *1000 / freq;
    }
    // ???????3
    double time4 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head4);
        iterativePlus_sum(a, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail4);
        time4 += (tail4 - head4) *1000 / freq;
    }
    // ??????
    cout << fixed << setprecision(6);
    cout << "?????? (n=" << n << ", trials=" << trials << "):\n";
    cout << "1. ???????????: " << time1 / trials << " ms\n";
    cout << "2. ???????????: " << time2 / trials << " ms\n";
    cout << "3. ?????????????: " << time3 / trials << " ms\n";
    cout << "4. ????????????: " << time4 / trials << " ms\n";

    // ???????4
    for (int t = 0; t < trials; ++t) {
        iterativePlus_sum(d, n);
    }


    // ??????
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;

    return 0;
}
