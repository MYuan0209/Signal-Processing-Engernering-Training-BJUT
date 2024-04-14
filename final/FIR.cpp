#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
#include"../func.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


const int table_pow_2[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};


//根据卷积公式编写卷积函数
//x[N]卷h[M],结果长度:L=M+N-1;
//由于:0<=k<=M+N-2,0<=n-k<=M-1,0<=k<=N-1
//所以:K>=0,K>=n-M+1;
//	 k<=N-1,K<=n; 
void convol(double x[], int N, double h[], int M, double y[], int L) {
	int k, kk;
	for(int n = 0; n < L; n++){
		for(k = (n >= M) ? n - M + 1 : 0, kk = (N > n) ? n : N - 1, y[n] = 0; k <= kk; k++) {
			y[n] += x[k] * h[n - k];
		}
	}
}


void convol_divided(double *x, int N, double *h, int M, double *y, int L, int divide_size) {
	int n_x;
	/*calculate the number of groups for x*/
	if (N % divide_size == 0) {
		n_x = N / divide_size;
	}
	else {
		n_x = (N + divide_size) / divide_size;
	}
	/*length after divided group*/
	int n, temp, i;
	for (i = 0; divide_size > table_pow_2[i]; ++i);
	temp = i;
	for (i = 0; M > table_pow_2[i]; ++i);
	n = max(temp, i);
	n = table_pow_2[n + 1];
	double *result = new double[n * n_x + n - 1];
	memset(result, 0, sizeof(result[0]) * (n * n_x + n - 1));
	/*update h*/
	COMPLEX *h_new = new COMPLEX[n];
	for (i = 0; i < M; ++i) {
		h_new[i].r = h[M - 1 - i];
		h_new[i].i = 0;
	}
	for (; i < n; ++i) {
		h_new[i].r = 0;
		h_new[i].i = 0;
	}
	fft(h_new, n, 1);
	for (i = 0; i < n_x; ++i) {
		int j;
		COMPLEX *xi = new COMPLEX[n];
		for (j = 0; j < divide_size; ++j) {
			xi[j].r = x[i * divide_size + j];
			xi[j].i = 0;
		}
		for (; j < n; ++j) {
			xi[j].r = 0;
			xi[j].i = 0;
		}
		fft(xi, n, 1);
		for (j = 0; j < n; ++j) {
			xi[j] *= h_new[j];
		}
		fft(xi, n, -1);
		for (j = 0; j < n; ++j) {
			result[i * divide_size + j] += xi[j].r;
		}
	}
	for (i = 0; i < L; ++i) {
		y[i] = result[i];
	}
}


void sampling (double *x, double (*f)(double, double), int N) {
	int T = 1;
	for (int i = 0; i < N; ++i) {
		x[i] = f(1.0 * i / N, T) + 0.05 * sin(2 * MY_PI * 450 * i / N);
	}
}


int main()
{
    double x[] = {1, 2, 3, 4, 5, 6, 7};
	double h[] = {3, 2, 1};
	int N = (sizeof(x) / sizeof(x[0])), M = (sizeof(h) / sizeof(h[0]));
	int L = N + M - 1;
	double *y1 = new double[L];
	double *y2 = new double[L];
	double *h1 = new double[M];
	for (int i = 0; i < M; ++i) {
	
		h1[i] = h[M - 1 - i];
	}
	cout << "convolution result:" << endl;
	convol(x, N, h1, M, y1, L);
	for (int i = 0; i < L; ++i) {
		cout << "y[" << i << "]=" << y1[i] << endl;
	}
	memset(y2, 0, sizeof(y2[0]) * L);
	convol_divided(x, N, h, M, y2, L, 2);
	cout << "divided convolution result" << endl;
	for (int i = 0; i < L; ++i) {
		cout << "y[" << i << "]=" << y2[i] << endl;
	}

	M = 50+1;
	int fs = 2000, f1 = 300, f2 = 600;
	double *b = new double[M];
	firDesgin(b, M, BANDSTOP, HAMMING, 1.0 * fs, 1.0 * f1, 1.0 * f2);
	cout << "print b:" << endl;
	for (int i = 0; i < M; ++i) {
		printf("%.6lf\t", b[i]);
		if ((i + 1) % 5 == 0) {
			cout << endl;
		}
	}

	window2((WCHAR*)"bandstop", 0, -80, fs / 2, 5, (char*)"hz", (char*)"db");
	xy2(BLUE);
	plotxy2(RED, 2, f, 20 * log10(firAbs(f, fs, b, M)));
	line2(f1, win2.y2, f1, win2.y1);
	line2(f2, win2.y2, f2, win2.y1);
	_getch();

	N = 2000;
	double *data = new double[N];
	delete []y1;
	delete []y2;
	sampling(data, rec, N);
	L = N + M - 1;
	y1 = new double[L];
	convol(data, N, b, M, y1, L);
	window2((WCHAR*)"FIR Digital Fitter", -0.5, -0.5, L, 2, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(WHITE, GREEN, data, N);
	_getch();
	plotgri2(WHITE, RED, y1, L);
	_getch();

	y2 = new double[L];
	convol_divided(data, N, b, M, y2, L, 50);
	plotgri2(WHITE, BLUE, y2, L);
	_getch();

	COMPLEX *X = new COMPLEX[N];
	COMPLEX *Y = new COMPLEX[L];
	COMPLEX *result_x = new COMPLEX[N];
	COMPLEX *result_y = new COMPLEX[L];
	for (int i = 0; i < N; ++i) {
		X[i].r = data[i];
		X[i].i = 0;
	}
	for (int i = 0; i < L; ++i) {
		Y[i].r = y1[i];
		Y[i].i = 0;
	}
	dft(result_x, X, N, 1);
	dft(result_y, Y, L, 1);
	for (int i = 0; i < N; ++i) {
		data[i] = abs(result_x[i]);
	}
	for (int i = 0; i < L; ++i) {
		y1[i] = abs(result_y[i]);
	}
	window2((WCHAR*)"frequency spectrum", 0, 0, L, 30, (char*)"k", (char*)"X[k]");
	xy2(BLACK);
	plot_frequency_spectrum(BLUE, data, N);
	_getch();
	frame2((char*)"k", (char*)"Y[k]");
	xy2(BLACK);
	plot_frequency_spectrum(BLUE, y1, L);
    return 0;
}
