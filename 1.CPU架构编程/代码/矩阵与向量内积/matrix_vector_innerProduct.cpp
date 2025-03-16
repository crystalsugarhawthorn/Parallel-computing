#include <iostream>
#include <cstdlib>
#include <Windows.h>
#include <iomanip>

// ����1�������ȱ��������治�Ѻã�
void method1(double** matrix, double* vector, double* result, int n) {
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += matrix[i][j] * vector[i]; // �����ȷ���
        }
        result[j] = sum;
    }
}

// ����2����ʽ�������ۼ�
void method2(double** b, double* a, double* sum, int n) {
    for (int i = 0; i < n; ++i) sum[i] = 0.0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            sum[i] += b[j][i] * a[j]; // �����ȷ���
        }
    }
}

// ����3��ѭ��չ��
void method3(double** b, double* a, double* sum, int n) {
    for (int i = 0; i < n; ++i) sum[i] = 0.0;
    for (int j = 0; j < n; ++j) {
    int i;
    // ��ѭ����ÿ�δ���4��Ԫ��
    for (i = 0; i <= n - 4; i += 4) {
         sum[i]     += b[j][i]     * a[j];
         sum[i + 1] += b[j][i + 1] * a[j];
         sum[i + 2] += b[j][i + 2] * a[j];
         sum[i + 3] += b[j][i + 3] * a[j];
    }
    // ����ʣ�಻��4����Ԫ��
    for (; i < n; ++i) {
         sum[i] += b[j][i] * a[j];
    }
}

}

int main() {
    const int n = 20000;      // ����ά�ȣ����Դ�ߴ�������
    const int trials = 1;   // ���Դ���

    // ��̬�����ڴ棨ģ�������ȴ洢��
    double** matrix = new double*[n]; // ����1�ľ��������ȴ洢��
    double** b = new double*[n];      // ����2�ľ��������ȴ洢��
    double* vector = new double[n];   // ��������
    double* a = new double[n];        // ������������
    double* result1 = new double[n];  // ����1���
    double* result2 = new double[n];  // ����2���
    double* result3 = new double[n];  // ����3���

    // �ֱ�Ϊ matrix �� b ����ÿһ�е��ڴ�
    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[n];
    }
    for (int i = 0; i < n; ++i) {
        b[i] = new double[n];
    }

    // ��ʼ������
    for (int i = 0; i < n; ++i) {
        vector[i] = static_cast<double>(rand()) / RAND_MAX;
        a[i] = vector[i];
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = static_cast<double>(rand()) / RAND_MAX;
            b[j][i] = matrix[i][j]; // �� b ��Ϊ matrix ��ת�ã�ģ�������ȣ�
        }
    }

    // �߾��ȼ�ʱ����
    long long head1, tail1, head2, tail2, freq, head3, tail3;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

    // ���Է���1
    double time1 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head1);
        method1(matrix, vector, result1, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);
        time1 += (tail1 - head1) *1000 / freq;
    }
    // ���Է���2
    double time2 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head2);
        method2(b, a, result2, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
        time2 += (tail2 - head2) *1000 / freq;
    }

    // ���Է���3
    double time3 = 0.0;
    for (int t = 0; t < trials; ++t) {
        QueryPerformanceCounter((LARGE_INTEGER*)&head3);
        method3(b, a, result3, n);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail3);
        time3 += (tail3 - head3) *1000 / freq;
    }

    // ������
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Method1 (�����ȷ���) ƽ��ʱ��: " << time1 / trials << " ����" << std::endl;
    std::cout << "Method2 (��ʽ���ۼ�) ƽ��ʱ��: " << time2 / trials << " ����" << std::endl;
    std::cout << "Method3 (ѭ���Ż�) ƽ��ʱ��: " << time3 / trials << " ����" << std::endl;

    // �ͷ��ڴ�
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
