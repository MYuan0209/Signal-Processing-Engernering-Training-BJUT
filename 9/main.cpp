#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
using namespace std;


extern WIN2 win2;
extern WIN3 win3;


double ahs1(double r,double i) {
	COMPLEX s(r, i), a;
	if (r == -1 && (i == 3 || i == -3)) {
		a = 0;
	}
	else {
		a = 10 * s * s / (s * s + 2 * s + 10);
	}
	return abs(a);
}


double phs1(double r,double i) {
	COMPLEX s = COMPLEX(r,i);
	COMPLEX a;
	a = 10 * s * s / (s * s + 2 * s + 10);
	return atan2(a.i, a.r);
}


double ahs2(double r,double i) {
	COMPLEX s(r, i), a;
	a = 200 / (s + 1000) / (s - 2000);
	return abs(a);
}


double phs2(double r,double i) {
	COMPLEX s = COMPLEX(r,i);
	COMPLEX a;
	a = 200 / (s + 1000) / (s - 2000);
	return atan2(a.i, a.r);
}


double ahs3(double r,double i) {
	COMPLEX s(r, i), a;
	a = 10 * (1 + s / 10) / (1 + s) / (1 + s / 100);
	return abs(a);
}


double phs3(double r,double i) {
	COMPLEX s = COMPLEX(r,i);
	COMPLEX a;
	a = 10 * (1 + s / 10) / (1 + s) / (1 + s / 100);
	return atan2(a.i, a.r);
}


//简单包装构造新的绘图工具：
#define _plotxyz3c(penColor,penWidth,x,y,fxy) \
{HPEN hPen=CreatePen(PS_SOLID, penWidth, penColor);\
HPEN hOldPen=(HPEN)SelectObject(win3.hdc, hPen);\
_plotxyz3(x,y,fxy);\
SelectObject(win2.hdc, hOldPen);\
DeleteObject(hPen);\
}


//简单包装构造新的绘图工具：
#define _plotxy2c(penColor,penWidth,x,fx) \
{HPEN hPen=CreatePen(PS_SOLID, penWidth, penColor);\
HPEN hOldPen=(HPEN)SelectObject(win3.hdc, hPen);\
_plotxy2(x,fx);\
SelectObject(win2.hdc, hOldPen);\
DeleteObject(hPen);\
}


int main() 
{
	window3((WCHAR*)"周期信号频谱图", -4, -2, -1, 8, 4, 16, (char*)"r", (char*)"i", (char*)"|h(s)|");
	xyz3(RED);
	_plotxyz3c(GREEN,1,r,i,ahs1(r,i));
	_getch();

	window2((WCHAR*)"2维幅频特性",-20, -10, 20, 10, (char*)"w", (char*)"|h(jw)|");
	_plotxy2c(RED,2,w,ahs1(0,w));
	_getch();

	window2((WCHAR*)"2维相频特性",-2,-3.15, 2, 3.15, (char*)"w", (char*)"phase");
	_plotxy2c(RED,2,w,phs1(0,w));
	_getch();

	window3((WCHAR*)"周期信号频谱图", -4, -2, -0.1, 8, 4, 0.1, (char*)"r", (char*)"i", (char*)"|h(s)|");
	xyz3(RED);
	_plotxyz3c(GREEN,1,r,i,ahs2(r,i));
	_getch();

	window2((WCHAR*)"2维幅频特性",-2, -1, 2, 2,   (char*)"w", (char*)"|h(jw)|");
	_plotxy2c(RED,2,w,ahs2(0,w));
	_getch();

	window2((WCHAR*)"2维相频特性",-2,-3.15,  2,3.15, (char*)"w", (char*)"phase");
	_plotxy2c(RED,2,w,phs2(0,w));
	_getch();

	window3((WCHAR*)"周期信号频谱图", -4, -2, -1, 8, 4, 16, (char*)"r", (char*)"i", (char*)"|h(s)|");
	xyz3(RED);
	_plotxyz3c(GREEN,1,r,i,ahs3(r,i));
	_getch();

	window2((WCHAR*)"2维幅频特性",-2, -1, 2, 2,   (char*)"w", (char*)"|h(jw)|");
	_plotxy2c(RED,2,w,ahs3(0,w));
	_getch();

	window2((WCHAR*)"2维相频特性",-2,-3.15,  2,3.15, (char*)"w", (char*)"phase");
	_plotxy2c(RED,2,w,phs3(0,w));
	_getch();
	return 0;
}