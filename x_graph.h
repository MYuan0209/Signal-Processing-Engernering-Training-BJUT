#ifndef X_GRAPH_H
#define X_GRAPH_H

#include<windows.h>
#include"x_math.h"

/*Global declaration of the GetConsoleWindow() function call.*/
extern "C" {WINBASEAPI HWND WINAPI GetConsoleWindow();}

#define WHITE RGB(255,255,255)
#define BLACK RGB(0,0,0)
#define RED RGB(255,0,0)
#define BLUE RGB(0,0,255)
#define GREEN RGB(0,0x80,0)
#define CYAN RGB(0,0xFF,0xFF)
#define MAGENTA RGB(0xFF,0,0xFF)
#define GRAY RGB(0x80, 0x80, 0x80)
#define LIGHTRED RGB(255,0xB6,0xC1)
#define LIGHTBLUE RGB(0xad,0xd8,0xe6)
#define LIGHTGREEN RGB(0x90,0xee,0x90)
#define LIGHTCYAN RGB(0xe0,0xFF,0xFF)
#define LIGHTMAGENTA RGB(0xFF,0x14,0x93)
#define LIGHTGRAY RGB(0xd3,0xd3,0xd3)
#define DARKGRAY RGB(0xa9,0xa9,0xa9)
#define YELLOW RGB(0xf0,0xe6,0x8c)

#define LRTB_SPACE 25 //Set the default top, bottom, left and right borders to white space.

const double ALPHA = MY_PI / 8;
const double BEITA = MY_PI / 6;


class WIN2 {
public:
	HWND hwnd;
	HDC hdc;
	RECT rect, rect2;
	char xstr[50], ystr[50];
	int xo, yo;
	double x1, y1, x2, y2;
	double xk, yk, dxPerPixel, dyPerPixel;
};


class WIN3 {
public:
	HWND hwnd;
	HDC hdc;
	RECT rect;
	char xstr[50], ystr[50], zstr[50];
	double x1, y1, z1, x2, y2, z2;
	double dxPerPixel, dyPerPixel, dzPerPixel;//dxPerPixel=1/xk,dyPerPixel=1/yk,dzPerPixel=1/zk.
	double xPixels, yPixels, zPixels, gridPixels;//xPixels is the length of the coordinate axis.
	double t[4][4];//t[0][] is row x, t[1][] is row y, t[2][] is row z, t[3] is constant row.
};


typedef union{
	struct {
		double x, y, z, c;
	};
	double A[4];
}POINT3,*PPOINT3;


typedef struct{
	int x, y;
}SCR2,*PSCR2;


class VERTEX{
public:
	double x;
	double y;
	double z;
	VERTEX();
	virtual ~VERTEX();
	VERTEX(double x1,double y1,double z1);
};


class EDGE{
public:
	int eNumber;
	int *eIndex;
    EDGE();
	virtual ~EDGE();
};


class box{
public:
	POINT3 v[8];
	double x0,y0,z0;
	double len;
	box();
	virtual ~box();
};


inline void showwin() {
	SendMessage(HWND_BROADCAST, WM_PAINT, 0, 0);
	_getch();
}


inline int ux2(double x, WIN2 win2) {
	return (int)(win2.xo+x*win2.xk);
}


inline int uy2(double y, WIN2 win2) {
	return (int)(win2.yo-y*win2.yk);
}


void draw_xyz(HPEN hpen1, HPEN hpen2, HPEN hpen3);

/*细去锯齿线，宽度为1。*/
void lineto1(HDC hdc,COLORREF bkcolor,COLORREF linecolor,int x,int y);

/*镶边线，将(x0,y0)到(x1,y1)直线的非(cx,cy)点所在边镶一个去锯齿边*/
void purfleline(HDC hdc,COLORREF bkcolor,COLORREF linecolor,int x0,int y0,int x1,int y1,int cx,int cy);
void purfleline(HDC hdc,COLORREF bkcolor,COLORREF linecolor,POINT a,POINT b,POINT ref);

/*箭头绘制函数，从点a到b的有向线段在b点绘一个箭头。*/
void arrowline(HDC hdc,int linewidth,COLORREF bkcolor,COLORREF linecolor,int ax,int ay,int bx,int by);
void arrowline(HDC hdc,int linewidth,COLORREF bkcolor,COLORREF linecolor,POINT a,POINT b);

/*画粗线（去锯齿）函数，width>=1.*/
void rectline(HDC hdc,int width,COLORREF bkcolor,COLORREF linecolor,double ax,double ay,double bx,double by);

void coord2(double x1,double y1,double x2,double y2);

void frame2(char xstr[],char ystr[]);
void frame2(char xstr[],char ystr[],int space);
void frame2();

void moveto2(double x,double y);
void lineto2(double x,double y);
void line2(double x1,double y1,double x2,double y2);

void plot2(double (*functionName)(double));
void plot2(double (*functionName)(double),double dt);

/*
inline void _plotxy2(double t, double (*f)(double), WIN2 win2) {
	t = win2.x1;
	moveto2(t, f(t));
	for (; t < win2.x2; t += win2.dxPerPixel) {
		lineto2(t, f(t));
	}
	t = win2.x2;
	lineto2(t, f(t));
}*/

#define _plotxy2(t,func) \
{double t;\
t=win2.x1;\
moveto2(t,func);\
for(;t<win2.x2;t=t+win2.dxPerPixel)lineto2(t,func);\
t=win2.x2;lineto2(t,func);}

/*The variable t does not need to be declared first*/
#define plotxy2(penColor,penWidth,t,func) \
{HPEN pen=CreatePen(PS_SOLID,penWidth,penColor),oldpen=(HPEN)SelectObject(win2.hdc,pen);\
_plotxy2(t,func);\
SelectObject(win2.hdc,oldpen);\
DeleteObject(pen);\
}


void moveto3(double x,double y,double z);
void moveto3(POINT3 a);

void lineto3(double x,double y,double z);
void lineto3(POINT3 a);

void line3(double x1,double y1,double z1,double x2,double y2,double z2);
void line3(POINT3 a,POINT3 b);

void mMul(double* A,double* B,int m,int n,int p,double* C);

POINT3 pt3(double x,double y,double z);

/*Translation coordinate*/
void shift3(double dx,double dy,double dz);
/*Translation of a set of points in the current coordinates.*/
void shift3(POINT3 p[],int n,double dx,double dy,double dz);

void rx3(POINT3 p[],int n,double a); //a is the angle of rotation.
void rx3(double a); //a is the rotation, the transformation of the X-axis.

void sxyz3(double sx,double sy,double sz);

void ry3(POINT3 p[],int n,double a); //a is the angle of rotation.
void ry3(double a); //a is the rotation, the transformation of the Y-axis.

void rz3(POINT3 p[],int n,double a); //a is the angle of rotation.
void rz3(double a); //a is the angle of rotation.

POINT getxy3();

void rv3(double w, double x1, double y1, double z1, double x2, double y2, double z2);
void rv3(POINT3 C[], int n, double w, double x1, double y1, double z1, double x2, double y2, double z2);

/*test success*/
void testRot();

void rv3();

void coord3(double x1,double y1,double z1,double x2,double y2,double z2);
void coord3(POINT3 a,POINT3 b);

/*bkcolor is not used, dotLinecolor is used and the background is white.*/
void frame3(COLORREF bkcolor,COLORREF dotLineColor,char *xstr,char *ystr,char *zstr);
void frame3(char *xstr,char *ystr,char *zstr);
void frame3();

void arrow3(COLORREF linecolor,double ax,double ay,double az,double bx,double by,double bz,int h,int d2);

void xyz3(COLORREF axiesColor);

void initgraph00(WCHAR title[],int x,int y,int width,int higth);

void initgraph(WCHAR title[],int x,int y,int width,int higth);


/*模仿函数方式，变量x,y需要先声明。*/
#define _plotxyz3(x,y,fxy) \
{double x,y,grid=win3.gridPixels*win3.dyPerPixel;\
for(y=win3.y1;y<=win3.y2;y=y+grid)\
{x=win3.x1,moveto3(x,y,fxy);\
for(x=win3.x1;x<=win3.x2;x=x+win3.dxPerPixel)\
lineto3(x,y,fxy);}\
y=win3.y2;\
{x=win3.x1,moveto3(x,y,fxy);\
for(x=win3.x1;x<=win3.x2;x=x+win3.dxPerPixel)\
lineto3(x,y,fxy);}\
grid=win3.gridPixels*win3.dxPerPixel;\
for(x=win3.x1;x<=win3.x2;x=x+grid)\
{y=win3.y1,moveto3(x,y,fxy);\
for(y=win3.y1;y<=win3.y2;y=y+win3.dyPerPixel)\
lineto3(x,y,fxy);}\
x=win3.x2;\
{y=win3.y1,moveto3(x,y,fxy);\
for(y=win3.y1;y<=win3.y2;y=y+win3.dyPerPixel)\
lineto3(x,y,fxy);}\
}

/*模仿函数方式，变量x,y需要先声明。*/
#define plotxyz3(penColor,penWidth,x,y,fxy) \
{HPEN hPen=CreatePen(PS_SOLID, penWidth, penColor);\
HPEN hOldPen=(HPEN)SelectObject(win3.hdc, hPen);\
_plotxyz3(x,y,fxy);\
SelectObject(win2.hdc, hOldPen);\
DeleteObject(hPen);\
}

void plot3(COLORREF color,double (*fxy)(double ,double));
void plot3(COLORREF gridcolor,COLORREF linecolor,COMPLEX p[],int M);

enum USEWIN{USE_CONSOLE,USE_WINDOW};

class winx{
public:
	HINSTANCE hinstance;
	HWND hwnd0;
	HDC hdc0;
	RECT rect0;
	winx(USEWIN win = USE_CONSOLE, char *title = (char*)"xxgc实验平台", int x0 = 0, int y0 = 0, int width = 800, int height = 600);
};

void xy2(COLORREF color);

void window2(WCHAR title[],double x1,double y1,double x2,double y2,char *xAxis,char *yAxis);

void clearviewport(COLORREF bkcolor);

void setvp3(int scrx,int scry,int width,int heigth,WCHAR title[],char *notice);

void window3(WCHAR title[],double x1,double y1,double z1,double x2,double y2,double z2,char *xAxis,char *yAxis,char *zAxis);

/*Clear the graphics of the 3 D user window.*/
void clear3();

void instKeyCtr3();

void instKeyCtr2();

/*Loop read keyboard input, hit ESC key, function end loop exit.
hit '+', '=', '.' or '>', the function increases the value in x[0] by the amount of (b-a)/n.
hit '_', '-', ',', or '<', the function reduces the number in x[0] by the amount of (b-a)/n.*/
BOOL keyCtr2(double x[1],double a,double b,int n);

BOOL keyCtr3(double x[1],double a,double b,int n);

/*Only clear the drawing area.*/
void clearvp2();

void instKeyCtr(double x[1],double a,double b,int N);

BOOL keyCtr(double x[1],double a,double b,int n);

/*draw 2D function through points*/
void plotgri2(COLORREF gridcolor, COLORREF linecolor, double *p, int N);

void plot_frequency_spectrum(COLORREF gridcolor, double *p, int N);

#endif