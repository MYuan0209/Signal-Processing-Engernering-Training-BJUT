#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
#include"../func.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


void rec(double* x, 
			double T, double N) {
	int i;
	for (i = 0; i < T; ++i) {
		x[i] = 1;
	}
	for (; i < N; ++i) {
		x[i] = 0;
	}
}


void fitter1(double* y, 
				double a1, double b0, double b1, double* x, int N) {
	y[0] += b0 * x[0];
	for (int i = 1; i < N; ++i) {
		y[i] += b0 * x[i] + b1 * x[i - 1] - a1 * y[i - 1];
	}
}


void fitter2(double* y, 
				double* a, double* b, double* x, int N) {
	y[0] += b[0] * x[0];
	y[1] += b[0] * x[1] + b[1] * x[0] - a[1] * y[0];
	for (int i = 2; i < N; ++i) {
		y[i] += b[0] * x[i] + b[1] * x[i - 1] + b[2] * x[i - 2] - a[1] * y[i - 1] - a[2] * y[i - 2];
	}
}


void fitter4(double* y, 
				double* a, double* b, double* x, int N) {
	y[0] += b[0] * x[0];
	y[1] += b[0] * x[1] + b[1] * x[0] - a[1] * y[0];
	y[2] += b[0] * x[2] + b[1] * x[1] + b[2] * x[0] - a[1] * y[1] - a[2] * y[0];
	y[3] += b[0] * x[3] + b[1] * x[2] + b[2] * x[1] + b[3] * x[0] - a[1] * y[2] - a[2] * y[1] - a[3] * y[0];
	for (int i = 4; i < N; ++i) {
		y[i] += b[0] * x[i] + b[1] * x[i - 1] + b[2] * x[i - 2] + b[3] * x[i - 3] + b[4] * x[i - 4] 
			- a[1] * y[i - 1] - a[2] * y[i - 2] - a[3] * y[i - 3] - a[4] * y[i - 4];
	}
}


int IIR_design(double H[][2][5], 
					int bandType, double db1, double db2, double fs, double f1, double f2, double f3, double f4) {
	int order, L;
	double b[10][2], c[2][3];

	order = btwOrder(bandType, db1, db2, fs, f1, f2, f3, f4);

	L = btwAf(b, order);

	btwC23(c, bandType, order, db1, fs, f1, f2, f3, f4);

	btwAf2Df(H, L, b, c);

	return L;
}


void IIR_draw_af(double H[][2][5], int L, int bandType, double db1, double db2, double fs, double f1, double f2, double f3, double f4) {
	window2((WCHAR*)"lowpass", -1., 5.,fs, -60., (char*)"hz", (char*)"db");
	xy2(BLUE);
	plotxy2(RED,2,f,btw20lgHz(f,fs,H,L));

	// draw reference line
	if(bandType==LOWPASS||bandType==HIGHPASS) {
		line2(0,-db1,win2.x2,-db1);
		line2(0,-db2,win2.x2,-db2);
		line2(f1,0,f1,win2.y1);
		line2(f2,0,f2,win2.y1);
	}
	else{
		line2(0,-db1,win2.x2,-db1);
		line2(0,-db2,win2.x2,-db2);
		line2(f1,0,f1,win2.y1);
		line2(f2,0,f2,win2.y1);
		line2(f3,0,f3,win2.y1);
		line2(f4,0,f4,win2.y1);
	}
	_getch();
}


int main()
{
	int N = 20, M=10;
	double *x, **y, **temp;
	x = new double[N];
	y = new double*[M];
	temp = new double*[M];
	for (int i = 0; i < M; ++i) {
		y[i] = new double[N];
		memset(y[i], 0, sizeof(y[i][0]) * N);
	}
	double a2[] = {1, 1, 0.25}, b2[] = {1, 2, 1};
	rec(x, N / 4, N);
	fitter1(y[1], 0.5, 1, 1, x, N);
	fitter1(y[2], 0.5, 1, 1, y[1], N);
	fitter2(y[0], a2, b2, x, 20);
	window2((WCHAR*)"IIR Digital Fitter", -0.5, -0.5, 20, 2.5, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(BLUE, GREEN, x, N);
	_getch();
	plotgri2(BLUE, RED, y[1], N);
	_getch();
	plotgri2(BLUE, CYAN, y[2], N);
	_getch();
	plotgri2(BLUE, MAGENTA, y[0], N);
	_getch();
	
	double a4[] = {1, 2, 1.5, 0.5, 0.0625}, b4[] = {1, 4, 6, 4, 1};
	memset(y[0], 0, sizeof(y[0][0]) * N);
	fitter2(y[3], a2, b2, y[2], N);
	fitter4(y[0], a4, b4, x, N);
	window2((WCHAR*)"IIR Digital Fitter", -0.5, -0.5, 20, 5, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(BLUE, GREEN, x, N);
	plotgri2(BLUE, CYAN, y[2], N);
	plotgri2(BLUE, LIGHTRED, y[3], N);
	_getch();
	plotgri2(BLUE, GRAY, y[0], N);
	_getch();

	int divide_size = 2;
	memset(y[0], 0, sizeof(y[0][0]) * N);
	double *x1 = new double[N];
	memset(x1, 0, sizeof(x1[0]) * N);
	for (int i = 0; i < N; i += divide_size) {
		for (int j = 0; j < divide_size; ++j) {
			x1[j] = x[i + j];
		}
		double *temp = new double[N - i];
		fitter4(temp, a4, b4, x1, N - i);
		for (int j = 0; j < N - i; ++j) {
			y[0][i + j] += temp[j];
		}
	}
	window2((WCHAR*)"IIR Digital Fitter", -0.5, -0.5, 20, 5, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(BLUE, GREEN, x, N);
	plotgri2(BLUE, LIGHTRED, y[3], N);
	_getch();
	plotgri2(BLUE, GRAY, y[0], N);
	_getch();

	delete []x;
	for (int i = 0; i < M; ++i) {
		delete []y[i];
	}
	delete []y;
	delete []x1;
	N = 1000;
	x = new double[N];
	x1 = new double[N];
	for (int j = 0; j < M; ++j) {
		temp[j] = new double[N];
	}
	double H[10][2][5];
	int L = IIR_design(H, 1, 3, 30, 1000, 100, 200, 300, 400);
	cout << "print H:" << endl;
	for (int i = 0; i < L; ++i) {
		cout << i << endl;
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 5; ++k) {
				printf("%.6lf\t", H[i][j][k]);
			}
			cout << endl;
		}
	}
 
	cout << "Normalization H:" << endl;
	for (int i = 0; i < L; ++i) {
		cout << i << endl;
		double temp = H[i][1][0];
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 5; ++k) {
				H[i][j][k] /= temp;
				printf("%.6lf\t", H[i][j][k]);
			}
			cout << endl;
		}
	}

	IIR_draw_af(H, L, 1, 3, 30, 1000, 100, 200, 300, 400);

	y = new double*[L];
	for (int i = 0; i < M; ++i) {
		y[i] = new double[N];
		memset(y[i], 0, sizeof(y[i][0]) * N);
	}
	rec(x, N / 2, N); // product a signal after sampling in a period
	fitter4(y[0], H[0][1], H[0][0], x, N);
	for (int i = 1; i < L; ++i) {
		fitter4(y[i], H[i][1], H[i][0], y[i - 1], N);
	}
	window2((WCHAR*)"IIR Digital Fitter", -0.5, -0.1, N, 1.5, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(WHITE, GREEN, x, N);
	_getch();
	plotgri2(WHITE, RED, y[L - 1], N);
	_getch();

	divide_size = 100;
	memset(y[0], 0, sizeof(y[0][0]) * N);
	memset(x1, 0, sizeof(x1[0]) * N);
	for (int i = 0; i < N; i += divide_size) {
		for (int j = 0; j < divide_size; ++j) {
			x1[j] = x[i + j];
		}
		for (int j = 0; j < L; ++j) {
			memset(temp[j], 0, sizeof(temp[j][0]) * (N - i));
		}
		fitter4(temp[0], H[0][1], H[0][0], x1, N - i);
		for (int j = 1; j < L; ++j) {
			fitter4(temp[j], H[j][1], H[j][0], temp[j - 1], N - i);
		}
		for (int j = 0; j < N - i; ++j) {
			y[0][i + j] += temp[L - 1][j];
		}
	}
	window2((WCHAR*)"IIR Digital Fitter", -0.5, -0.1, N, 1.5, (char*)"n", (char*)"y[n]");
	xy2(BLACK);
	plotgri2(WHITE, GREEN, x, N);
	plotgri2(WHITE, RED, y[L - 1], N);
	_getch();
	plotgri2(WHITE, CYAN, y[0], N);
	_getch();

	COMPLEX *Y = new COMPLEX[N];
	COMPLEX *X = new COMPLEX[N];
	COMPLEX *result_x = new COMPLEX[N];
	COMPLEX *result_y = new COMPLEX[N];
	for (int i = 0; i < N; ++i) {
		X[i] = COMPLEX(x[i], 0);
		Y[i] = COMPLEX(y[L - 1][i], 0);
	}
	dft(result_x, X, N, 1);
	dft(result_y, Y, N, 1);
	for (int i = 0; i < N; ++i) {
		x[i] = sqrt(result_x[i].r * result_x[i].r + result_x[i].i * result_x[i].i);
		y[L - 1][i] = sqrt(result_y[i].r * result_y[i].r + result_y[i].i * result_y[i].i);
	}
	window2((WCHAR*)"IIR Digital Fitter", 0, 0, N, 25, (char*)"k", (char*)"X[k]");
	xy2(BLACK);
	plot_frequency_spectrum(BLUE, x, N);
	_getch();
	frame2((char*)"k", (char*)"Y[k]");
	xy2(BLACK);
	plot_frequency_spectrum(BLUE, y[L - 1], N);
	return 0;
}