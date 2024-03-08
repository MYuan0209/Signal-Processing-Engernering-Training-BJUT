#include<bits/stdc++.h>
#include<windows.h>
#include<conio.h>
#include"../x_math.h"
#include"../x_graph.h"
#include"func.h"
using namespace std;


int main()
{
	extern HWND hWndFather;
	extern WIN2 win2;
	extern WIN3 win3;

	double t;

  	window2((WCHAR*)"Function graphic display",-8, 6, 8, -6, (char*)"t",(char*)"f(t)");

	HPEN pen1=CreatePen(PS_SOLID,3,RED), oldpen=(HPEN)SelectObject(win2.hdc,pen1);
	HPEN pen2=CreatePen(PS_SOLID,3,BLUE);
	//plot the cosine function
	xy2(RED);
	_plotxy2(t, cos(t));
	_getch();
	//plot the tangent fuction
	frame2(win2.xstr,win2.ystr);	
	xy2(RGB(200,120,0));
	SelectObject(win2.hdc, pen1);
	_plotxy2(t, tan(t));
	_getch();
	//plot the exponential function
	frame2(win2.xstr,win2.ystr);	
	xy2(RGB(200,120,0));
	SelectObject(win2.hdc, pen1);
	_plotxy2(t, exp(t));
	_getch();
	//plot the ramp function
	frame2(win2.xstr,win2.ystr);	
	xy2(RGB(200,120,0));
	SelectObject(win2.hdc, pen1);
	_plotxy2(t, ramp(t));
	_getch();
	//plot the sawtooth wave
	frame2(win2.xstr,win2.ystr);	
	xy2(RGB(200,120,0));
	SelectObject(win2.hdc, pen1);
	_plotxy2(t, sawtoothwave(t, 4));
	_getch();
	//plot the equation x * x + y * y = 4
	frame2(win2.xstr, win2.ystr);
	xy2(RGB(200, 120, 0));
	moveto2(2.0, 0.0);
	for (int i = 0; i < 360000; ++i) {
		SelectObject(win2.hdc, pen1);
		lineto2(2 * cos(1.0 * i / 180000 * MY_PI), 2 * sin(1.0 * i / 180000 * MY_PI));
		moveto2(2 * cos(1.0 * i / 180000 * MY_PI), 2 * sin(1.0 * i / 180000 * MY_PI));
	}
	_getch();
	/*plot the trangle*/
	frame2(win2.xstr, win2.ystr);
	xy2(RGB(200, 120, 0));
	moveto2(0, 0);
	SelectObject(win2.hdc, pen1);
	lineto2(3, 3);
	lineto2(3, 0);
	lineto2(0, 0);
	_getch();
	//plot the sector
	frame2(win2.xstr, win2.ystr);
	xy2(RGB(200, 120, 0));
	moveto2(2.0, 0.0);
	for (int i = 0; i < 120; ++i) {
		SelectObject(win2.hdc, pen1);
		lineto2(2 * cos(1.0 * i / 180 * MY_PI), 2 * sin(1.0 * i / 180 * MY_PI));
		moveto2(2 * cos(1.0 * i / 180 * MY_PI), 2 * sin(1.0 * i / 180 * MY_PI));
	}
	lineto2(0.0, 0.0);
	line2(0.0, 0.0, 2.0, 0.0);
	_getch();

	window3((WCHAR*)"三维函数显示图",-3,-3,-10,3,3,10,(char*)"x",(char*)"y",(char*)"z");
	/*plot the sinc function*/
	xyz3(RED);
	SelectObject(win3.hdc,pen2);
	_plotxyz3(t, w, sinc(t, w, 1000));
	_getch();

	SelectObject(win2.hdc,oldpen);
	DeleteObject(pen1);
	DeleteObject(pen2);
  	return 0;
}
