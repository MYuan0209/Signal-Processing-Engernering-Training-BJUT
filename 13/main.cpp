#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


void rec(double* x, double T, double N) {
	int i;
	for (i = 0; i < T; ++i) {
		x[i] = 1;
	}
	for (; i < N; ++i) {
		x[i] = 0;
	}
}


void tri(double* x, double T, double N) {
	int i;
	for (i = 0; i <= T; ++i) {
		x[i] = 0.2 * i;
	}
	for (; i <= 2 * T; ++i) {
		x[i] = 2 - 0.2 * i;
	}
	for (; i < N; ++i) {
		x[i] = 0;
	}
}


void my_sin(double* x, double T, double N) {
	int i;
	for (i = 0; i < T; ++i) {
		x[i] = sin(2 * MY_PI / T * i);
	}
	for (; i < N; ++i) {
		x[i] = 0;
	}
}


void delta(double* x, double N) {
	x[0] = 1;
	for (int i = 1; i < N; ++i) {
		x[i] = 0;
	}
}


/*void plotgri2(COLORREF gridcolor, COLORREF linecolor, double *p, int N) {
	HPEN pen1 = CreatePen(PS_SOLID,1,gridcolor), oldpen=(HPEN)SelectObject(win3.hdc,pen1);
	HPEN pen2 = CreatePen(PS_SOLID,1,linecolor);
	for(int i = 0; i < N; i++) {
		line2(i, 0.0, i, p[i]);
	}

	SelectObject(win3.hdc, pen2);
	moveto2(0, p[0]);
	for(int i = 0; i < N; i++) {
		lineto2(i, p[i]);
	}

	SelectObject(win2.hdc, oldpen);
	DeleteObject(pen1);
	DeleteObject(pen2);
}*/


void fitter1(double* y, double a0, double b0, double b1, double* x) {
	const int N = 20;
	y[0] = b0 * x[0];
	for (int i = 1; i < 20; ++i) {
		y[i] = b0 * x[i] + b1 * x[i - 1] + a0 * y[i - 1];
	}
}


void fitter2(double* y, double* a, double* b, double* x) {
	const int N = 20;
	y[0] = b[0] * x[0];
	y[1] = b[0] * x[1] + b[1] * x[0] + a[0] * y[0];
	for (int i = 2; i < N; ++i) {
		y[i] = b[0] * x[i] + b[1] * x[i - 1] + b[2] * x[i - 2] + a[0] * y[i - 1] + a[1] * y[i - 2];
	}
}


int main() 
{
	int N = 20;
	double x[20], y[20], y2[20], temp[20];
	double a[5], b[5];

	window2((WCHAR*)"Digital Fitter",-2, -2, 20, 5, (char*)"n", (char*)"y(n)");

	rec(x, 11, 20);
	fitter1(y, 0.4, 1, 0.2, x);
	xy2(RED);
	plotgri2(BLUE, GREEN, x, N);
	plotgri2(BLUE, RED, y, N);
	_getch();

	tri(x, 5, 20);
	fitter1(y, 0.4, 1, 0.2, x);
	frame2((char*)"n", (char*)"y(n)");
	xy2(RED);
	plotgri2(BLUE, GREEN, x, N);
	plotgri2(BLUE, RED, y, N);
	_getch();

	my_sin(x, 16, 20);
	fitter1(y, 0.4, 1, 0.2, x);
	frame2((char*)"n", (char*)"y(n)");
	xy2(RED);
	plotgri2(BLUE, GREEN, x, N);
	plotgri2(BLUE, RED, y, N);
	_getch();

	delta(x, 20);
	a[0] = 1.5;
	a[1] = -0.5;
	b[0] = 0;
	b[1] = 1;
	b[2] = 0;
	fitter1(temp, 0, 1, 0.5, x);
	fitter1(y2, 1, 0, 1, temp);
	fitter2(y, a, b, x);
	frame2((char*)"n", (char*)"y(n)");
	xy2(RED);
	plotgri2(BLUE, GREEN, x, N);
	//plotgri2(BLUE, YELLOW, temp, N);
	//plotgri2(BLUE, RED, y2, N);
	_getch();
	plotgri2(BLUE, CYAN, y, N);
	_getch();

	return 0;
}
