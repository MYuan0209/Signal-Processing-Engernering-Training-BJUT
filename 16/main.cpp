#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


double squarewave(double t, double T, double h) {
	double a;
	if (T <= 0 || h <= 0) {
		printf("squareWave() input error\n");
		exit(1);
	}
	a = fmod(t, T);
	if (a > 0) {
		return (a <= h) ? 1.0 : 0.0;
	}
	else {
		return (a > (-T + h)) ? 0.0 : 1.0;
	} 
}


void convol(double x[], int N, double h[], int M, double y[], int L) {
	int n, k, kk;
	for(n = 0; n < L; n++){
		for(k=(n>=M)?n-M+1:0,kk=(N>n)?n:N-1,y[n]=0;k<=kk;k++) {
			y[n]+=x[k]*h[n-k];
		}
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


int main()
{
	//设计低通FIR
	//滤波器长度N=21;
	int N = 21;
	double b[21], fl = 100, fs = 1000;
	firDesgin(b, 21, LOWPASS, RECTANG, fs, fl);
	double x[400], y[420];
	for (int i = 0; i < 400; ++i) {
		x[i] = squarewave(1.0 * i / 400, 0.01, 0.005);
	}
	conv(y, x, 400, b, 21);

	window2((WCHAR*)"example", -1, -2, 400, 2, (char*)"x", (char*)"y");
	xy2(BLUE);
	plotxy2(GREEN, 2, t, squarewave(1.0 * t / 400, 400, 0.005));
	// plotgri2(GREEN, GREEN, x, 400);
	plotgri2(RED, RED, y, 420);
	_getch();
	/*
	double b[21], fl = 100, fs = 1000;
	firDesgin(b, 21, LOWPASS, RECTANG, fs, fl);

	window2((WCHAR*)"低通幅频特性", -fs / 2, 2, fs / 2, -2, (char*)"hz", (char*)"|H|");
	xy2(BLUE);
	plotxy2(GREEN, 2, f, firAbs(f, fs, b, N));
	_getch();

	window2((WCHAR*)"低通相频特性", -fs / 2, 5, fs / 2, -5, (char*)"hz", (char*)"FI");
	xy2(BLUE);
	plotxy2(RED, 3, f, firPhase(f, fs, b, N));
	_getch();*/
/*
	{
		double b[21],fl=100,fh=400,fs=1000;
		firDesgin(b,21,BANDPASS,RECTANG,fs,fl,fh);
		window2((WCHAR*)"带通幅频特性",-fs/2,2,fs/2,-2,(char*)"hz",(char*)"|H|");
		xy2(BLUE);
		plotxy2(GREEN,2,f,firAbs(f,fs,b,N));
		_getch();

		window2((WCHAR*)"带通相频特性",-fs/2,5,fs/2,-5,(char*)"hz",(char*)"FI");
		xy2(BLUE);
		plotxy2(RED,3,f,firPhase(f,fs,b,N));
		_getch();
	}*/
	return 0;
}
