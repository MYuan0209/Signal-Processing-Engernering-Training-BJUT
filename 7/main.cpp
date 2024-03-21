#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


double winfunc(double t,double wide)//实际是个方波
{
	return (0<=t&&t<wide)?1:0;
}


double triwave(double t, double T, double wide) {
	if (t >= 0 && t < wide) {
		return t / wide;
	}
	else if (t >= wide && t < T) {
		return 1 - (t - wide) / (T - wide);
	}
	else {
		return 0;
	}
}


double aphsine(double t, double T, double h) {
	double result = sin(2 * MY_PI / T * t) - h + 1;
	if (t >= 0 && t < T) {
		return result > 0 ? result : 0;
	}
	return 0;
}


//构建绘制周期谱的工具函数
void plotgri3(COLORREF gridcolor, COLORREF linecolor, COMPLEX *p, int N) {
	HPEN pen1 = CreatePen(PS_SOLID,1,gridcolor), oldpen=(HPEN)SelectObject(win3.hdc,pen1);
	HPEN pen2 = CreatePen(PS_SOLID,1,linecolor);
	for(int i = 0; i < N; i++) {
		line3(i, 0, 0, i, p[i].r, p[i].i);
	}

	SelectObject(win3.hdc, pen2);
	moveto3(0, p[0].r, p[0].i);
	for(int i = 0; i < N; i++) {
		lineto3(i, p[i].r, p[i].i);
	}

	SelectObject(win2.hdc, oldpen);
	DeleteObject(pen1);
	DeleteObject(pen2);
}


int main() 
{
	int i, N = 20;
	COMPLEX p[20]; 
	
	window3((WCHAR*)"周期信号频谱图",-1,-0.2,-0.2,21,0.2,0.2, (char*)"N", (char*)"r", (char*)"i");
	
    xyz3(BLUE);
	double wide, T = 4;
	for(wide=0, instKeyCtr(&wide, 0, 4,20); keyCtr(&wide, 0, 4,20);) {
		for(i = 0; i < N; i++) {
			/*
			_iTrapz(p[i].r,0.0,T,10,winfunc(t,wide)*cos(MY_PI*2.0/T*t*i),t);
			p[i].r=p[i].r*2/T;
			_iTrapz(p[i].i,0.0,T,10,winfunc(t,wide)*sin(MY_PI*2.0/T*t*i),t);
			p[i].i=-p[i].i*2/T;
			*//*
			_iTrapz(p[i].r, 0.0, T, 10, triwave(t, T, wide) * cos(MY_PI * 2.0 / T * t * i), t);
			p[i].r = p[i].r * 2 / T;
			_iTrapz(p[i].i, 0.0, T, 10, triwave(t, T, wide) * sin(MY_PI * 2.0 / T * t * i), t);
			p[i].i = -p[i].i * 2 / T;
			*/
			_iTrapz(p[i].r, 0.0, T, 10, aphsine(t, T, wide) * cos(MY_PI * 2.0 / T * t * i), t);
			p[i].r = p[i].r * 2 / T;
			_iTrapz(p[i].i, 0.0, T, 10, aphsine(t, T, wide) * sin(MY_PI * 2.0 / T * t * i), t);
			p[i].i = -p[i].i * 2 / T;
		}

		p[0].r=p[0].r/2;//处理直流项
		frame3((char*)"N", (char*)"r", (char*)"i");
		instKeyCtr(&wide,0,4,20);
		xyz3(BLUE);
		plotgri3(BLUE,RED,p,N);
	}
	return 0;
}