#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


double sin_abs(double t) {
	if (t >= 0 && t < MY_PI) {
		return sin(t);
	}
	else {
		return 0;
	}
}


double triwave(double t) {
	if (t >= 0 && t < 2) {
		return t;
	}
	else if (t >= 2 && t < 4) {
		return 4 - t;
	}
	else {
		return 0;
	}
}


double squarewave(double t) {
	if (t >= 0 && t < 2) {
		return 1;
	}
	else {
		return 0;
	}
}


int main()
{
	extern WIN2 win2;
	extern WIN3 win3;
	double  a[3][20], b[3][20], c[3][20], t, T1 = MY_PI, T2 = 4;
	int i, N = 5;

	printf("\n积分运算需要时间,请等待!\n");

	for (int i = 0; i < N; ++i) {
		_iTrap(a[0][i], 0.0, T1, sin_abs(t) * cos(MY_PI * 2.0 / T1 * t * i), t);
		a[0][i] = a[0][i] * 2 / T1;
		_iTrap(b[0][i], 0.0, T1, sin_abs(t) * sin(MY_PI * 2.0 / T1 * t * i), t);
		b[0][i] = b[0][i] * 2 / T1;
		_iTrap(c[0][i], 0.0, T1, sin_abs(t) * wal(t / T1, i), t);
		c[0][i] = c[0][i] / T1;
	}
	a[0][0] = a[0][0] / 2;

	for (int i = 0; i < N; ++i) {
		_iTrap(a[1][i], 0.0, T2, triwave(t) * cos(MY_PI * 2.0 / T2 * t * i), t);
		a[1][i] = a[1][i] * 2 / T2;
		_iTrap(b[1][i], 0.0, T2, triwave(t) * sin(MY_PI * 2.0 / T2 * t * i), t);
		b[1][i] = b[1][i] * 2 / T2;
		_iTrap(c[1][i], 0.0, T2, triwave(t) * wal(t / T2, i), t);
		c[1][i] = c[1][i] / T2;
	}
	a[1][0] = a[1][0] / 2;

	for (int i = 0; i < N; ++i) {
		_iTrap(a[2][i], 0.0, T2, squarewave(t) * cos(MY_PI * 2.0 / T2 * t * i), t);
		a[2][i] = a[2][i] * 2 / T2;
		_iTrap(b[2][i], 0.0, T2, squarewave(t) * sin(MY_PI * 2.0 / T2 * t * i), t);
		b[2][i] = b[2][i] * 2 / T2;
		_iTrap(c[2][i], 0.0, T2, squarewave(t) * wal(t / T2, i), t);
		c[2][i] = c[2][i] / T2;
	}
	a[2][0] = a[2][0] / 2;

	/*
	for(i=0;i<N;i++){_fcosCoef(a[i],i,T,jieti(t),t);
	_fsinCoef(b[i],i,T,jieti(t),t);
	_WalCoef(c[i],i,T,jieti(t),t)}
	*/  
	printf("ok!积分运算结束，敲任意键盘继续\n");

	window2((WCHAR*)"训练3",-15, 2, 15, -2, (char*)"t", (char*)"f");
	xy2(BLUE);
	
	HPEN pen1=CreatePen(PS_SOLID,3,RED), oldpen=(HPEN)SelectObject(win2.hdc,pen1);
	HPEN pen2=CreatePen(PS_SOLID,1,GREEN);

	_plotxy2(t, FouSer(t, T1, a[0], b[0], N));
	_getch();

	SelectObject(win2.hdc, pen2);
	_plotxy2(t, sin_abs(t));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	SelectObject(win2.hdc,pen1);
	_plotxy2(t, WalSer(t, T1, c[0], N));
	_getch();

	SelectObject(win2.hdc,pen2);
	_plotxy2(t, sin_abs(t));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	SelectObject(win2.hdc,pen1);
	_plotxy2(t, sin_abs(t) * sin_abs(t) - FouSer(t, T1, a[0], b[0], N) * FouSer(t, T1, a[0], b[0], N));
	SelectObject(win2.hdc,pen2);
	_plotxy2(t, sin_abs(t) * sin_abs(t) - WalSer(t, T1, c[0], N) * WalSer(t, T1, c[0], N));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	_plotxy2(t, FouSer(t, T2, a[1], b[1], N));
	_getch();

	SelectObject(win2.hdc, pen2);
	_plotxy2(t, triwave(t));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	SelectObject(win2.hdc,pen1);
	_plotxy2(t, WalSer(t, T2, c[1], N));
	_getch();

	SelectObject(win2.hdc,pen2);
	_plotxy2(t, triwave(t));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	SelectObject(win2.hdc,pen1);
	_plotxy2(t, triwave(t) * triwave(t) - FouSer(t, T2, a[1], b[1], N) * FouSer(t, T2, a[1], b[1], N));
	SelectObject(win2.hdc,pen2);
	_plotxy2(t, triwave(t) * triwave(t) - WalSer(t, T2, c[1], N) * WalSer(t, T2, c[1], N));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	_plotxy2(t, FouSer(t, T2, a[2], b[2], N));
	_getch();

	SelectObject(win2.hdc, pen2);
	_plotxy2(t, squarewave(t));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	SelectObject(win2.hdc,pen1);
	_plotxy2(t, WalSer(t, T2, c[2], N));
	_getch();

	SelectObject(win2.hdc,pen2);
	_plotxy2(t, squarewave(t));
	_getch();

	frame2(win2.xstr, win2.ystr);
	xy2(RED);
	SelectObject(win2.hdc,pen1);
	_plotxy2(t, squarewave(t) * squarewave(t) - FouSer(t, T2, a[2], b[2], N) * FouSer(t, T2, a[2], b[2], N));
	SelectObject(win2.hdc,pen2);
	_plotxy2(t, squarewave(t) * squarewave(t) - WalSer(t, T2, c[2], N) * WalSer(t, T2, c[2], N));
	_getch();

	SelectObject(win2.hdc, oldpen);//恢复原来画笔.
	DeleteObject(pen1);//删除新建画笔.
	DeleteObject(pen2);//删除新建画笔.

	return 0;
}
