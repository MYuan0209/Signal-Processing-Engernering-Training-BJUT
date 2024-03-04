#include<bits/stdc++.h>
#include <windows.h>  //This lib MUST be included 
#include <conio.h>
#include <tchar.h>
#include"x_math.h"
#include"x_graph.h"
using namespace std;


HWND hWndFather;
WIN2 win2;
WIN3 win3;


VERTEX::VERTEX() {

}


VERTEX::~VERTEX() {

}


VERTEX::VERTEX(double x1, double y1, double z1) {
	x=x1;
	y=y1;
	z=z1;
}


void draw_xyz(HPEN hpen1, HPEN hpen2, HPEN hpen3) {
	extern WIN3 win3;
	POINT ca;
	SelectObject(win3.hdc,hpen1);
	line3(0,0,0,40,0,0);
	ca=getxy3();
	TextOutW(win3.hdc,ca.x,ca.y,(LPCWSTR)"x",1);
	showwin();
	SelectObject(win3.hdc,hpen2);
	line3(0,0,0,0,40,0);
	ca=getxy3();
	TextOutW(win3.hdc,ca.x,ca.y,(LPCWSTR)"y",1);
	showwin();
	SelectObject(win3.hdc,hpen3);
	line3(0,0,0,0,0,40);
	ca=getxy3();
	TextOutW(win3.hdc,ca.x,ca.y,(LPCWSTR)"z",1);
	showwin();
}


void lineto1(HDC hdc,COLORREF bkcolor,COLORREF linecolor,int x,int y) {
	POINT currentPosition;
	MoveToEx(hdc, x, y, &currentPosition);
	int x0 = currentPosition.x, y0 = currentPosition.y;

	if (x == x0) {
		if (y0 > y) {
			swap(y, y0);
		}
		for(; y0 < y; y0++) {
			SetPixelV(hdc, x0, y0, linecolor);
		}
	}
	else if (y == y0) {
		if (x0 > x) {
			swap(x0, x);
		}
		for(; x0 < x; x0++) {
			SetPixelV(hdc, x0, y0, linecolor);
		}
	}
	else {
		int rf = GetRValue(linecolor), gf = GetGValue(linecolor), bf = GetBValue(linecolor);
		int re = GetRValue(bkcolor) - rf, ge = GetGValue(bkcolor) - gf, be = GetBValue(bkcolor) - bf;
		double k = (double)(y - y0) / (double)(x - x0), k2 = 1.0 / k, t; //Calculated slope increment
		if(k > 1.0 || k < -1.0) { //y principal direction
			if(y0 > y) {
				swap(x0, x);
				swap(y0, y);
			}		
			for(t = (double)x0; y0 <= y; y0++, t+=k2) {
				double e1 = t - (int)t, e2 = 1 - e1;
				SetPixelV(hdc, (int)t, y0, RGB(round(re*e1)+rf,round(ge*e1)+gf,round(be*e1)+bf));
				SetPixelV(hdc, (int)t+1, y0, RGB(round(re*e2)+rf,round(ge*e2)+gf,round(be*e2)+bf));
			}
		}
		if(-1.0 <= k && k <= 1.0){//x principal direction
			if(x0 > x) {
				swap(x0, x);
				swap(y0, y);
			}		
			for(t=y0;x0<=x;x0++,t+=k) {
				double e1=t-(int)t,e2=1-e1;
				SetPixelV(hdc,x0,(int)t,RGB(round(re*e1)+rf,round(ge*e1)+gf,round(be*e1)+bf));
				SetPixelV(hdc,x0,(int)t+1,RGB(round(re*e2)+rf,round(ge*e2)+gf,round(be*e2)+bf));
			}	
		}
	}
}


void purfleline(HDC hdc,COLORREF bkcolor,COLORREF linecolor,int x0,int y0,int x1,int y1,int cx,int cy) {
	auto min = [] (int a, int b) {
		return a < b ? a : b;
	};

	double t0, e0;
	/*if k2 equals to +-INF, then k=0.*/
	double k2 = (double)(y1 - y0) / (double)(x1 - x0), k = -1.0 / k2;
	int rf = GetRValue(linecolor), gf = GetGValue(linecolor), bf = GetBValue(linecolor);
	int re = GetRValue(bkcolor) - rf, ge = GetGValue(bkcolor) - gf, be = GetBValue(bkcolor)-bf;

	if(-1.0>k2||k2>1.0) { //y principal direction
		if(cx-x0-k*(cy-y0)>0) {
			if(y0>y1) {
				swap(x0, x1);
				swap(y0, y1);
			}
			for(t0=(double)x0;y0<=y1;y0++,t0-=k) {
				e0=t0-(int)t0;
				COLORREF c=GetPixel(hdc,(int)t0,y0);
				int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
				SetPixelV(hdc,(int)t0,y0,RGB(r,g,b));
				SetPixelV(hdc,(int)t0+1,y0,linecolor);
			}
		}
		else {
			if(y0>y1) {
				swap(x0, x1);
				swap(y0, y1);
			}
			for(t0=(double)x0;y0<=y1;y0++,t0-=k) {
				e0=1-(t0-(int)t0);
				SetPixelV(hdc,(int)t0,y0,linecolor);
				COLORREF c=GetPixel(hdc,(int)t0+1,y0);
				int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
				SetPixelV(hdc,(int)t0+1,y0,RGB(r,g,b));
			}
		}
	}
	if(-1.0<=k2&&k2<=1.0) { //x principal direction
		if(cy-y0-k2*(cx-x0)>0) {
			if(x0>x1) {
				swap(x0, x1);
				swap(y0, y1);
			}
			for(t0=(double)y0;x0<=x1;x0++,t0+=k2) {
				e0=t0-(int)t0;
				COLORREF c=GetPixel(hdc,x0,(int)t0);
				int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
				SetPixelV(hdc,x0,(int)t0,RGB(r,g,b));
				SetPixelV(hdc,x0,(int)t0+1,linecolor);
			}//end for
		}
		else {
			if(x0>x1) {
				swap(x0, x1);
				swap(y0, y1);
			}
			for(t0=(double)y0;x0<=x1;x0++,t0+=k2) {
				e0=1-t0+(int)t0;
				SetPixelV(hdc,x0,(int)t0,linecolor);
				COLORREF c=GetPixel(hdc,x0,(int)t0+1);
				int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
				SetPixelV(hdc,x0,(int)t0+1,RGB(r,g,b));
			}
		}	
	}
}//end purfleline


void purfleline(HDC hdc,COLORREF bkcolor,COLORREF linecolor,POINT a,POINT b,POINT ref) {
	purfleline(hdc,bkcolor,linecolor,a.x,a.y,b.x,b.y,ref.x,ref.y);
}


void arrowline(HDC hdc,int linewidth,COLORREF bkcolor,COLORREF linecolor,int ax,int ay,int bx,int by) {
	int arrowbottomwidth=3*linewidth,arrowhigth=4*linewidth;
	double k,k2;
	/*v[0],v[1],v[2],v[6] points form a rectangle counterclockwise. 
	v[3],v[4],v[5] form the arrow counterclockwise, and v[4] is the tip of the arrow*/
	POINT v[7];
	if(ay==by) {
		v[0].x=ax,v[0].y=ay+linewidth/2;
		v[1].x=ax,v[1].y=ay-linewidth/2;
		v[6].x=bx,v[6].y=v[0].y;
		v[2].x=bx,v[2].y=v[1].y;
		v[5].x=bx,v[5].y=by+arrowbottomwidth/2;
		v[3].x=bx,v[3].y=by-arrowbottomwidth/2;
		v[4].x=bx+((bx>ax)?arrowhigth:-arrowhigth);v[4].y=by;
	}
	else if(ax==bx) {
		v[0].x=ax+linewidth/2,v[0].y=ay;
		v[1].x=ax-linewidth/2,v[1].y=ay;
		v[6].x=bx+linewidth/2,v[6].y=by;
		v[2].x=bx-linewidth/2,v[2].y=by;
		v[5].x=bx+arrowbottomwidth/2,v[5].y=by,
		v[3].x=bx-arrowbottomwidth/2,v[3].y=by;
		v[4].x=bx,v[4].y=by+((by>ay)?arrowhigth:-arrowhigth);
	}
	else {
		k2=(double)(ax-bx)/(double)(ay-by);
		//(x-x0)sqrt(1+k2*k2)=+-wide/2.0
		double st=linewidth/(2.0*sqrt(1+k2*k2));
		v[0].x=round(ax+st);v[0].y=round(ay-k2*st);
		v[1].x=round(ax-st);v[1].y=round(ay+k2*st);
		v[6].x=round(bx+st);v[6].y=round(by-k2*st);
		v[2].x=round(bx-st);v[2].y=round(by+k2*st);
		st=arrowbottomwidth/(2.0*sqrt(1+k2*k2));
		v[5].x=round(bx+st);v[5].y=round(by-k2*st);
		v[3].x=round(bx-st);v[3].y=round(by+k2*st);
		k=-1/k2;
		st=arrowhigth/(sqrt(1+k*k));
		v[4].x=round(bx+((bx>ax)?st:-st));
		v[4].y=round(by-k*((bx>ax)?st:-st));
	}
	
	HBRUSH brush=CreateSolidBrush(linecolor);
	HBRUSH oldbrush=(HBRUSH)SelectObject(hdc,brush);
	HPEN pen=CreatePen(PS_SOLID,1,linecolor);
	HPEN oldpen=(HPEN)SelectObject(hdc,pen);

	Polygon(hdc,v,7);

	int rf=GetRValue(linecolor),gf=GetGValue(linecolor),bf=GetBValue(linecolor);
	int re=GetRValue(bkcolor)-rf,ge=GetGValue(bkcolor)-gf,be=GetBValue(bkcolor)-bf;
	purfleline(hdc,bkcolor,linecolor,v[0],v[6],v[1]);
	purfleline(hdc,bkcolor,linecolor,v[1],v[2],v[0]);
	//���ƶ̱ߣ�/////////////////////////////////////
	purfleline(hdc,bkcolor,linecolor,v[0],v[1],v[6]);
	//����3�㵽6����߶Σ��б�ӱ߷���ʱ��Ҫ�ߵ���
	purfleline(hdc,bkcolor,linecolor,v[6],v[5],v[4]);
	//���Ƶ�4����5���߱�,�б�ӱ߷���ʱ��Ҫ�ߵ���
	purfleline(hdc,bkcolor,linecolor,v[2],v[3],v[4]);
	//���Ƶ�6����7���߱�
	purfleline(hdc,bkcolor,linecolor,v[3],v[4],v[5]);
	//���Ƶ�5����7���߱�
	purfleline(hdc,bkcolor,linecolor,v[4],v[5],v[3]);

	SelectObject(hdc,oldpen);
	DeleteObject(pen);
	SelectObject(hdc,oldbrush);
	DeleteObject(brush);
}


void arrowline(HDC hdc,int linewidth,COLORREF bkcolor,COLORREF linecolor,POINT a,POINT b) {
	arrowline(hdc,linewidth,bkcolor,linecolor,a.x,a.y,b.x,b.y);
}


void rectline(HDC hdc,int width,COLORREF bkcolor,COLORREF linecolor,double ax,double ay,double bx,double by) {
	auto min = [] (int a, int b) {
		return a < b ? a : b;
	};

	double k,k2;
	double x1,y1,x2,y2,x3,y3,x4,y4;
	POINT v[4], vt;

	if(fabs(ay-by)<1.0e-6) {
		x1=ax,y1=ay+width/2.0;
		x2=ax,y2=ay-width/2.0;
		x3=bx,y3=y1;
		x4=bx,y4=y2;
	}
	else {
		k2=-(ax-bx)/(ay-by);
		x1=ax+width/(2.0*sqrt(1+k2*k2));y1=ay+k2*width/(2*sqrt(1+k2*k2));
		x2=ax-width/(2.0*sqrt(1+k2*k2));y2=ay-k2*width/(2*sqrt(1+k2*k2));
		x3=bx+width/(2.0*sqrt(1+k2*k2));y3=by+k2*width/(2*sqrt(1+k2*k2));
		x4=bx-width/(2.0*sqrt(1+k2*k2));y4=by-k2*width/(2*sqrt(1+k2*k2));
	}

	/*11,12,--22,21 points form a rectangle. 
	Line x11x21(v[0]v[3]) is parallel to line x12x22(v[1]v[2]). 
	Line x11x12(v[0]v[1]) is parallel to line 21-22(v[3]v[2].*/
	v[0].x=round(x1);v[0].y=round(y1);
	v[1].x=round(x2);v[1].y=round(y2);
	v[2].x=round(x4);v[2].y=round(y4);
	v[3].x=round(x3);v[3].y=round(y3);

	HBRUSH brush=CreateSolidBrush(linecolor);
	HBRUSH oldbrush=(HBRUSH)SelectObject(hdc,brush);
	Polygon(hdc,v,4);
	SelectObject(hdc,oldbrush);
	DeleteObject(brush);

	int rf=GetRValue(linecolor),gf=GetGValue(linecolor),bf=GetBValue(linecolor);
	int re=GetRValue(bkcolor)-rf,ge=GetGValue(bkcolor)-gf,be=GetBValue(bkcolor)-bf;

	double t0,t1,e0,e1;

	if(v[0].x!=v[3].x&&v[0].y!=v[3].y){
		k=(double)(v[3].y-v[0].y)/(double)(v[3].x-v[0].x),k2=1.0/k;
		if(k>1.0||k<-1.0) { //y principal direction
			if(v[0].y > v[3].y) {
				swap(v[0], v[3]);
				swap(v[1], v[2]);
			}
			if(v[0].x > v[1].x) {
				swap(v[0], v[1]);
				swap(v[3], v[2]);
			}
			for(t0=(double)v[0].x,t1=(double)v[1].x;v[0].y<=v[3].y;t0+=k2,t1+=k2,v[0].y++,v[1].y++) {
				e0=t0-(int)t0,e1=1-e0;
				COLORREF c=GetPixel(hdc,(int)t0,v[0].y);
				int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
				SetPixelV(hdc,(int)t0,v[0].y,RGB(r,g,b));
				SetPixelV(hdc,(int)t0+1,v[0].y,linecolor);

				SetPixelV(hdc,(int)t1,v[1].y,linecolor);

				c=GetPixel(hdc,(int)t1+1,v[0].y);
				r=min(round(re*e1)+rf,GetRValue(c)),g=min(round(ge*e1)+gf,GetGValue(c)),b=min(round(be*e1)+bf,GetBValue(c));
				SetPixelV(hdc,(int)t1+1,v[1].y,RGB(r,g,b));
			}
		}
		if(-1<=k&&k<=1.0) { //x principal direction
			if(v[0].x>v[3].x) {
				swap(v[0], v[3]);
				swap(v[1], v[2]);
			}
			if(v[0].y>v[1].y) {
				swap(v[0], v[1]);
				swap(v[3], v[2]);
			}
			for(t0=(double)v[0].y,t1=(double)v[1].y;v[0].x<=v[3].x;t0+=k ,t1+=k ,v[0].x++,v[1].x++) {
				e0=t0-(int)t0,e1=1-e0;
				COLORREF c=GetPixel(hdc,v[0].x,(int)t0);
				int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
				SetPixelV(hdc,v[0].x,(int)t0,RGB(r,g,b));
				SetPixelV(hdc,v[0].x,(int)t0+1,linecolor);
				SetPixelV(hdc,v[1].x,(int)t1,linecolor);
				c=GetPixel(hdc,v[1].x,(int)t1+1);
				r=min(round(re*e1)+rf,GetRValue(c)),g=min(round(ge*e1)+gf,GetGValue(c)),b=min(round(be*e1)+bf,GetBValue(c));
				SetPixelV(hdc,v[1].x,(int)t1+1,RGB(r,g,b));
			}
		}

		/*draw short edge*/
		v[0].x=round(x1);v[0].y=round(y1);
		v[1].x=round(x2);v[1].y=round(y2);
		v[2].x=round(x4);v[2].y=round(y4);
		v[3].x=round(x3);v[3].y=round(y3);
		vt=v[0],v[0]=v[3],v[3]=v[2],v[2]=v[1],v[1]=vt;
		k=(double)(v[3].y-v[0].y)/(double)(v[3].x-v[0].x),k2=1.0/k;
		if(k>1.0||k<-1.0) { //y principal direction
			if(v[0].y>v[3].y){vt=v[0],v[0]=v[3],v[3]=vt,vt=v[1],v[1]=v[2],v[2]=vt;}
			if(v[0].x>v[1].x){vt=v[0],v[0]=v[1],v[1]=vt,vt=v[3],v[3]=v[2],v[2]=vt;}
			for(t0=(double)v[0].x,t1=(double)v[1].x;v[0].y<=v[3].y;t0+=k2,t1+=k2,v[0].y++,v[1].y++)
			{e0=t0-(int)t0,e1=1-e0;
			COLORREF c=GetPixel(hdc,(int)t0,v[0].y);
			int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
			SetPixelV(hdc,(int)t0,v[0].y,RGB(r,g,b));
			SetPixelV(hdc,(int)t0+1,v[0].y,linecolor);
			//SetPixelV(hdc,(int)t0+2,v[0].y,RGB(255,0,0));
			//SetPixelV(hdc,(int)t1-2,v[1].y,RGB(0,255,0));		
			c=GetPixel(hdc,(int)t1+1,v[0].y);
			r=min(round(re*e1)+rf,GetRValue(c)),g=min(round(ge*e1)+gf,GetGValue(c)),b=min(round(be*e1)+bf,GetBValue(c));
			SetPixelV(hdc,(int)t1,v[1].y,linecolor);
			SetPixelV(hdc,(int)t1+1,v[1].y,RGB(r,g,b));
			}
		}
		if(-1<=k&&k<=1.0) { //x principal direction
			if(v[0].x>v[3].x){vt=v[0],v[0]=v[3],v[3]=vt,vt=v[1],v[1]=v[2],v[2]=vt;}
			if(v[0].y>v[1].y){vt=v[0],v[0]=v[1],v[1]=vt,vt=v[3],v[3]=v[2],v[2]=vt;}
			for(t0=(double)v[0].y,t1=(double)v[1].y;v[0].x<=v[3].x;t0+=k ,t1+=k ,v[0].x++,v[1].x++)
			{e0=t0-(int)t0,e1=1-e0;
			COLORREF c=GetPixel(hdc,v[0].x,(int)t0);
			int r=min(round(re*e0)+rf,GetRValue(c)),g=min(round(ge*e0)+gf,GetGValue(c)),b=min(round(be*e0)+bf,GetBValue(c));
			SetPixelV(hdc,v[0].x,(int)t0,RGB(r,g,b));
			SetPixelV(hdc,v[0].x,(int)t0+1,linecolor);

			//SetPixelV(hdc,v[1].x,(int)t1-2,RGB(255,0,0));		
			c=GetPixel(hdc,v[1].x,(int)t1+1);
			r=min(round(re*e1)+rf,GetRValue(c)),g=min(round(ge*e1)+gf,GetGValue(c)),b=min(round(be*e1)+bf,GetBValue(c));
			SetPixelV(hdc,v[1].x,(int)t1,linecolor);
			SetPixelV(hdc,v[1].x,(int)t1+1,RGB(r,g,b));
			}
		}
	}
	/*This method of stroke is very poor, so image filtering should be considered.*/
}


void coord2(double x1,double y1,double x2,double y2) {
	if (x1>x2){win2.x2=x1;win2.x1=x2;}
	else if(x1==x2){printf("横坐标设置错误！");_getch();exit(0);}
	else {win2.x1=x1;win2.x2=x2;}

	if(y1>y2){win2.y2=y1;win2.y1=y2;}
	else if(y1==y2){printf("纵坐标设置错误！");_getch();exit(0);}
	else {win2.y1=y1;win2.y2=y2;}

	win2.rect2.left=win2.rect.left+LRTB_SPACE;
	win2.rect2.right=win2.rect.right-LRTB_SPACE;
	win2.rect2.top=win2.rect.top+LRTB_SPACE;
	win2.rect2.bottom=win2.rect.bottom-LRTB_SPACE;

	win2.xk=(double)(win2.rect2.right-win2.rect2.left)/(win2.x2-win2.x1);
	win2.yk=(double)(win2.rect2.bottom-win2.rect2.top)/(win2.y2-win2.y1); 
	win2.dxPerPixel=1.0/win2.xk;
	win2.dyPerPixel=1.0/win2.yk;

	win2.xo=(int)(-win2.x1*win2.xk)+win2.rect2.left;
	win2.yo=(int)(win2.y2*win2.yk+win2.rect2.top);
}


void frame2(char xstr[],char ystr[]) {
	const int DGRID = 20; //Number of pixels per grid

	SelectObject(win2.hdc,(HPEN)GetStockObject(BLACK_PEN));
	SelectObject(win2.hdc,(HBRUSH)GetStockObject(WHITE_BRUSH));//LTGRAY_BRUSH));
	Rectangle(win2.hdc,win2.rect.left,win2.rect.top,win2.rect.right,win2.rect.bottom);//ʹmemdcָ���λͼ��ˢ��huiɫ�� 

	WCHAR str[100]={0};//WCHAR�ǿ��ַ�����UNICODE��ʽ����
	swprintf_s(str,100,L"(%g,%g)",win2.x1,win2.y1);//100�ǻ���str��С
	TextOutW(win2.hdc,win2.rect2.left,win2.rect2.bottom+6,(LPCWSTR)str,wcslen(str));
	swprintf_s(str,100,L"(%g,%g)",win2.x2,win2.y2);
	TextOutW(win2.hdc,win2.rect2.right-8*wcslen(str)+14,win2.rect2.top-18,(LPCWSTR)str,wcslen(str));
	swprintf_s(str,100,L"%S:%g/grid",ystr,DGRID/(win2.yk));
	TextOutW(win2.hdc,win2.rect2.left,win2.rect2.top-18,(LPCWSTR)str,wcslen(str));
	swprintf_s(str,100,L"%S:%g/grid",xstr,DGRID/(win2.xk));//S��д����xstr��ANSI�ַ���
	TextOutW(win2.hdc,win2.rect2.right-wcslen(str)*8,win2.rect2.bottom+6,(LPCWSTR)str,wcslen(str));

	int i;
	for(i=win2.xo;i<win2.rect2.right;i=i+DGRID) {
		if(i>win2.rect2.left) {
			MoveToEx(win2.hdc,i,win2.rect2.bottom,NULL);
			LineTo(win2.hdc,i,win2.rect2.bottom+5);
		}
	}

	for(i=win2.xo;i>win2.rect2.left;i=i-DGRID) {
		if(i<win2.rect2.right) {
			MoveToEx(win2.hdc,i,win2.rect2.bottom,NULL);
			LineTo(win2.hdc,i,win2.rect2.bottom+5);
		}
	}

	for(i=win2.yo;i<win2.rect2.bottom;i=i+DGRID) {
		if(i>win2.rect2.top) {
			MoveToEx(win2.hdc,win2.rect2.left,i,NULL);
			LineTo(win2.hdc,win2.rect2.left-5,i);
		}
	}

	for(i=win2.yo;i>win2.rect2.top;i=i-DGRID) {
		if(i<win2.rect2.bottom) {
			MoveToEx(win2.hdc,win2.rect2.left,i,NULL);
			LineTo(win2.hdc,win2.rect2.left-5,i);
		}
	}

	if(win2.xo>=win2.rect2.left&&win2.xo<=win2.rect2.right) {
		MoveToEx(win2.hdc,win2.xo,win2.rect2.bottom,NULL);
		LineTo(win2.hdc,win2.xo-2,win2.rect2.bottom+5);
		LineTo(win2.hdc,win2.xo+2,win2.rect2.bottom+5);
		LineTo(win2.hdc,win2.xo,win2.rect2.bottom);
	}

	if(win2.yo>=win2.rect2.top&&win2.yo<=win2.rect2.bottom) {
		MoveToEx(win2.hdc,win2.rect2.left,win2.yo,NULL);
		LineTo(win2.hdc,win2.rect2.left-5,win2.yo-2);
		LineTo(win2.hdc,win2.rect2.left-5,win2.yo+2);
		LineTo(win2.hdc,win2.rect2.left,win2.yo);
	}
	
	Rectangle(win2.hdc,win2.rect2.left-1,win2.rect2.top-1,win2.rect2.right+1,win2.rect2.bottom+1);

	sprintf_s(win2.xstr,50,xstr);
	sprintf_s(win2.ystr,50,ystr);
}


void frame2(char xstr[],char ystr[],int space) {
	win2.rect2.left=win2.rect.left+space;
	win2.rect2.right=win2.rect.right-space;
	win2.rect2.top=win2.rect.top+space;
	win2.rect2.bottom=win2.rect.bottom-space;

	win2.xk=(double)(win2.rect2.right-win2.rect2.left)/(win2.x2-win2.x1);
	win2.yk=(double)(win2.rect2.bottom-win2.rect2.top)/(win2.y2-win2.y1); 
	win2.dxPerPixel=1.0/win2.xk;
	win2.dyPerPixel=1.0/win2.yk;
	win2.xo=(int)(-win2.x1*win2.xk)+win2.rect2.left;
	win2.yo=(int)(win2.y2*win2.yk+win2.rect2.top);

	frame2(xstr,ystr);
}


void frame2() {
	frame2(win2.xstr,win2.ystr);
}


void moveto2(double x,double y) {
	MoveToEx(win2.hdc,ux2(x, win2),uy2(y, win2),NULL);
}


void lineto2(double x,double y) {
	LineTo(win2.hdc,ux2(x, win2),uy2(y, win2));
}


void line2(double x1,double y1,double x2,double y2) {
	moveto2(x1,y1);
	lineto2(x2,y2);
}


void plot2(double (*functionName)(double)) {
	double x;
	x=win2.x1;
	moveto2(x,(double)(*functionName)(x));
	for(;x<win2.x2;x=x+win2.dxPerPixel)lineto2(x,(double)(*functionName)(x));
}


void plot2(double (*functionName)(double),double dt) {
	double x;
	x=win2.x1;
	moveto2(x,(double)(*functionName)(x));
	for(;x<win2.x2;x=x+dt)lineto2(x,(double)(*functionName)(x));
	x=win2.x2;lineto2(x,(double)(*functionName)(x));
}


void moveto3(double x,double y,double z) {
	int sx,sy;
	sx=(int)(x*win3.t[0][0]+y*win3.t[1][0]+z*win3.t[2][0]+win3.t[3][0]+0.5);
	sy=(int)(-x*win3.t[0][1]-y*win3.t[1][1]-z*win3.t[2][1]-win3.t[3][1]+0.5);
	MoveToEx(win3.hdc,sx,sy,NULL);
}


void moveto3(POINT3 a) {
	moveto3(a.x,a.y,a.z);
}


void lineto3(double x,double y,double z) {
	int sx,sy;
	sx=(int)(x*win3.t[0][0]+y*win3.t[1][0]+z*win3.t[2][0]+win3.t[3][0]+0.5);
	sy=(int)(-x*win3.t[0][1]-y*win3.t[1][1]-z*win3.t[2][1]-win3.t[3][1]+0.5);
	LineTo(win3.hdc,sx,sy);
}


void lineto3(POINT3 a) {
	lineto3(a.x,a.y,a.z);
}


void line3(double x1,double y1,double z1,double x2,double y2,double z2) {
	moveto3(x1,y1,z1);
	lineto3(x2,y2,z2);
}


void line3(POINT3 a,POINT3 b) {
	moveto3(a);
	lineto3(b);
}


void mMul(double* A,double* B,int m,int n,int p,double* C)
{
	int i,j,k;
	for(i=0;i<m;i++) {
		for(j=0;j<p;j++) {
			*(C+i*p+j)=0.0;
			for(k=0;k<n;k++) {
				(*(C+i*p+j))+=(*(A+i*n+k))*(*(B+k*p+j));
			}
		}
	}
}


POINT3 pt3(double x,double y,double z) {
	POINT3 a;
	a.x=x;
	a.y=y;
	a.z=z;
	a.A[3]=1;
	return a; 
}


void shift3(double dx,double dy,double dz) {
	double RX[4][4]={0};
	RX[0][0]=1;
	RX[1][1]=1;
	RX[2][2]=1;
	RX[3][3]=1;
	RX[3][0]=dx;
	RX[3][1]=dy;
	RX[3][2]=dz;

	double temp[4][4];
	mMul((double*)RX,(double*)win3.t,4,4,4,(double*)temp);
	int i,j;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			win3.t[i][j]=temp[i][j];
		}
	}
}


void shift3(POINT3 p[],int n,double dx,double dy,double dz) {
	int i;
	for(i=0;i<n;i++){
		p[i].x=p[i].x-dx;
		p[i].y=p[i].y-dy;
		p[i].z=p[i].z-dz;
	}
}


void rx3(POINT3 p[],int n,double a) {
	double RX[4][4]={0};
	RX[0][0]=1;
	RX[1][1]=cos(a);
	RX[2][2]=RX[1][1];
	RX[1][2]=sin(a);
	RX[2][1]=-RX[1][2];
	RX[3][3]=1;

	POINT3 temp;
	int i;
	for(i=0;i<n;i++){
		p[i].A[3]=1;
		mMul((double*)p[i].A,(double*)RX,1,4,4,(double*)temp.A);
		p[i].x=temp.x;p[i].y=temp.y;p[i].z=temp.z;
	}

}


void rx3(double a) {
	double RX[4][4]={0};
	RX[0][0]=1;
	RX[1][1]=cos(a);
	RX[2][2]=RX[1][1];
	RX[1][2]=sin(a);
	RX[2][1]=-RX[1][2];
	RX[3][3]=1;

	double temp[4][4];
	mMul((double*)RX,(double*)win3.t,4,4,4,(double*)temp);
	int i,j;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			win3.t[i][j]=temp[i][j];
		}
	}
}


void sxyz3(double sx,double sy,double sz) {
	double ST[4][4]={0};
	ST[0][0]=sx;
	ST[1][1]=sy;
	ST[2][2]=sz;
	ST[3][3]=1;

	double temp[4][4];
	mMul((double*)ST,(double*)win3.t,4,4,4,(double*)temp);
	int i,j;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			win3.t[i][j]=temp[i][j];
		}
	}
}


void ry3(POINT3 p[],int n,double a) {
	double RX[4][4]={0};
	RX[0][0]=cos(a);
	RX[1][1]=1;
	RX[2][2]=RX[0][0];
	RX[2][0]=sin(a);
	RX[0][2]=-RX[2][0];
	RX[3][3]=1;

	POINT3 temp;
	int i;
	for(i=0;i<n;i++){
		p[i].A[3]=1;
		mMul((double*)p[i].A,(double*)RX,1,4,4,(double*)temp.A);
		p[i].x=temp.x;
		p[i].y=temp.y;
		p[i].z=temp.z;
	}
}


void ry3(double a) {
	double RX[4][4]={0};
	RX[0][0]=cos(a);
	RX[1][1]=1;
	RX[2][2]=RX[0][0];
	RX[2][0]=sin(a);
	RX[0][2]=-RX[2][0];
	RX[3][3]=1;

	double temp[4][4];
	mMul((double*)RX,(double*)win3.t,4,4,4,(double*)temp);
	int i,j;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			win3.t[i][j]=temp[i][j];
		}
	}
}


void rz3(POINT3 p[],int n,double a) {
	double RX[4][4]={0};
	RX[0][0]=cos(a);
	RX[0][1]=sin(a);
	RX[1][0]=-RX[0][1];
	RX[1][1]=RX[0][0];
	RX[2][2]=1;
	RX[3][3]=1;

	POINT3 temp;
	int i;
	for(i=0;i<n;i++){
		p[i].A[3]=1;
		mMul((double*)p[i].A,(double*)RX,1,4,4,(double*)temp.A);
		p[i].x=temp.x;
		p[i].y=temp.y;
		p[i].z=temp.z;
	}
}


void rz3(double a) {
	double RX[4][4]={0};
	RX[0][0]=cos(a);
	RX[0][1]=sin(a);
	RX[1][0]=-RX[0][1];
	RX[1][1]=RX[0][0];
	RX[2][2]=1;
	RX[3][3]=1;

	double temp[4][4];
	mMul((double*)RX,(double*)win3.t,4,4,4,(double*)temp);
	int i,j;
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			win3.t[i][j]=temp[i][j];
		}
	}
}


POINT getxy3() {
	POINT currentPosition;
	MoveToEx(win3.hdc, 0, 0, &currentPosition);
	MoveToEx(win3.hdc, currentPosition.x, currentPosition.y, NULL);
	return currentPosition;
}


void rv3(double w,double x1,double y1,double z1, double x2,double y2,double z2) {
	double x=x2-x1, y=y2-y1, z=z2-z1;
	double wy, wz;

	/*The Angle of rotation of the vector about the y axis to the xoy plane*/
	wy=atan2(z,x);

	/*The Angle of rotation of the vector about the z axis to the x axis in the xoy plane*/
	wz=atan2(y,sqrt(x*x+z*z));

	shift3(x1,y1,z1);
	ry3(wy);
	rz3(-wz);
	rx3(w);
	rz3(wz);
	ry3(-wy);
	shift3(-x1,-y1,-z1);
}


void rv3(POINT3 C[], int n, double w, double x1, double y1, double z1, double x2, double y2, double z2) {
	double x=x2-x1,y=y2-y1,z=z2-z1;
	double wy,wz;

	/*The Angle of rotation of the vector about the y axis to the xoy plane*/
	wy=atan2(z,x);

	/*The Angle of rotation of the vector about the z axis to the x axis in the xoy plane*/
	wz=atan2(y,sqrt(x*x+z*z));

	shift3(C,n,x1,y1,z1);
	ry3(C,n,wy);
	rz3(C,n,-wz);
	rx3(C,n,w);
	rz3(C,n,wz);
	ry3(C,n,-wy);
	shift3(C,n,-x1,-y1,-z1);
}


void testRot() {
	HPEN hpen1 = CreatePen(PS_SOLID, 1, RGB(255, 0, 0));
	HPEN hpen2 = CreatePen(PS_SOLID, 1, RGB(0, 255, 0));
	HPEN hpen3 = CreatePen(PS_SOLID, 1, RGB(0, 0, 255));
	HPEN hpen4 = CreatePen(PS_SOLID, 1, RGB(200, 0, 255));
	HPEN hpen5 = CreatePen(PS_SOLID, 1, RGB(0, 200, 255));
	HPEN hpen6 = CreatePen(PS_SOLID, 1, RGB(200, 200, 0));

	rx3(BEITA);
	ry3(BEITA);
	draw_xyz(hpen1, hpen2, hpen3);

	POINT3 C[1]={20,100,0,1};

	POINT3 A={20,100,0,1};
	double x1=0,y1=0,z1=0,x2=0,y2=0,z2=100;
	int i;

	line3(x1,y1,z1,x2,y2,z2);

	moveto3(A);
	for(i=0;i<41;i++){
	rv3(C,1,3.1415926/20,x1,y1,z1,x2,y2,z2);
	lineto3(C[0]);}

	showwin();
	for(i=0;i<21;i++){
	rv3(3.1415926/10,0,100,0,200,0,0);
	draw_xyz(hpen1, hpen2, hpen3);}
	showwin();
}


void rv3() {
	HPEN hpen1=CreatePen(PS_SOLID,3,RGB(255,0,0));
	HPEN hpen2=CreatePen(PS_SOLID,3,RGB(0,255,0));
	HPEN hpen3=CreatePen(PS_SOLID,3,RGB(0,0,255));
	HPEN hpen4=CreatePen(PS_SOLID,3,RGB(200,0,255));
	HPEN hpen5=CreatePen(PS_SOLID,3,RGB(0,200,255));
	HPEN hpen6=CreatePen(PS_SOLID,3,RGB(200,200,0));

	rx3(BEITA);
	ry3(BEITA);
	shift3(-100,-100,-100);

	draw_xyz(hpen1, hpen2, hpen3);

	POINT3 c[3]={20,20,20,1,100,100,100,1,100,100,200,1};
	POINT3 org={0,0,0};
	double x=-100,y=100,z=100,wy,wz;

	auto drawline = [] (WIN3 win3, HPEN hpen, POINT3 c[]) {
		SelectObject(win3.hdc, hpen);
		line3(c[0], c[1]);
		lineto3(c[2]);
		showwin();
	};

	drawline(win3, hpen5, c);
	shift3(c,3,20,20,20);
	drawline(win3, hpen6, c);

	x=c[1].x;
	y=c[1].y;
	z=c[1].z;

	/*The Angle of rotation of the vector about the y axis to the xoy plane*/
	wy=atan2(z,x);

	/*The Angle of rotation of the vector about the z axis to the x axis in the xoy plane*/
	wz=atan2(y,sqrt(x*x+z*z));

	ry3(c,3,wy);
	drawline(win3, hpen5, c);

	rz3(c,3,-wz);
	printf("\ndd=%f,%f,%f;\n",c[1].x,c[1].y,c[1].z);
	drawline(win3, hpen5, c);

	/*After the vector has been rotated to the X-axis, 
	find the new coordinates of the point c at the Angle of rotation w about the vector.*/
	rx3(c,3,BEITA);

	rz3(c,3,wz);
	drawline(win3, hpen4, c);

	ry3(c,3,-wy);
	drawline(win3, hpen4, c);

	shift3(c,3,-20,-20,-20);
	drawline(win3, hpen5, c);

	double x1=20,y1=30,z1=40,x2=200,y2=-200,z2=80,xp=-20,yp=-100,zp=50; 
	SelectObject(win3.hdc,hpen4);line3(x1,y1,z1,x2,y2,z2);

	SelectObject(win3.hdc,hpen5);lineto3(xp,yp,zp);

	shift3(x1,y1,z1);
}


void coord3(double x1,double y1,double z1,double x2,double y2,double z2) {
	if (x1>x2) {
		win3.x2=x1;
		win3.x1=x2;
	}
	else if(x1==x2) {
		printf("x坐标设置错误!");
		_getch();
		exit(1);
	}
	else {
		win3.x1=x1;
		win3.x2=x2;
	}

	if(y1>y2) {
		win3.y2=y1;
		win3.y1=y2;
	}
	else if(y1==y2) {
		printf("y坐标设置错误!");
		_getch();
		exit(1);
	}
	else {
		win3.y1=y1;
		win3.y2=y2;
	}

	if(z1>z2) {
		win3.z2=z1;
		win3.z1=z2;
	}
	else if(z1==z2) {
		printf("z坐标设置错误!");
		_getch();
		exit(1);
	}
	else {
		win3.z1=z1;
		win3.z2=z2;
	}

	GetClientRect(win3.hwnd, &win3.rect);

	double w, h;
	w=win3.rect.right-win3.rect.left;
	h=-win3.rect.top+win3.rect.bottom;

	/*Since the following coordinate system first rotates 30 degrees around the X-axis 
	and then -30 degrees around the Y-axis, only the Y-axis projection is parallel to the window vertical. 
	The maximum side length of a coordinate equilateral stereo box is only half the vertical height of the window:*/
	if(h>w)h=w;
	win3.xPixels=h/2.0;

	win3.zPixels=win3.yPixels=win3.xPixels=(double)(10.0*((int)win3.xPixels/10));
	win3.dxPerPixel=(win3.x2-win3.x1)/win3.xPixels;
	win3.dyPerPixel=(win3.y2-win3.y1)/win3.yPixels;
	win3.dzPerPixel=(win3.z2-win3.z1)/win3.zPixels;

	/*Set the coordinate system to screen physical coordinates, 
	where the origin of the coordinates is in the upper left corner of the drawing client area*/
	win3.t[0][0]=1;win3.t[0][1]=0;win3.t[0][2]=0;win3.t[0][3]=0;
	win3.t[1][0]=0;win3.t[1][1]=1;win3.t[1][2]=0;win3.t[1][3]=0;
	win3.t[2][0]=0;win3.t[2][1]=0;win3.t[2][2]=0;win3.t[2][3]=0;
	win3.t[3][0]=0;win3.t[3][1]=0;win3.t[3][2]=0;win3.t[3][3]=1;

	/*Translate the origin of the physical coordinates to the center of the drawing client*/
	shift3((win3.rect.right+win3.rect.left)/2.0,-(win3.rect.top+win3.rect.bottom)/2.0,0);

	rx3(MY_PI/6);
	ry3(-MY_PI/6);
	rx3(-MY_PI/2);

	/*Transform user coordinates proportionally to physical coordinates, 
	so that window coordinates become user coordinate system*/
	sxyz3(1.0/win3.dxPerPixel,1.0/win3.dyPerPixel,1.0/win3.dzPerPixel);

	/*Under user coordinates, move and adjust the origin position of coordinates*/
	shift3(-(win3.x1+win3.x2)/2,-(win3.y1+win3.y2)/2,-(win3.z1+win3.z2)/2);

	win3.gridPixels=5;
}


void coord3(POINT3 a,POINT3 b) {
	coord3(a.x,a.y,a.z,b.x,b.y,b.z);
}


void frame3(COLORREF bkcolor,COLORREF dotLineColor,char *xstr,char *ystr,char *zstr) {
	SelectObject(win3.hdc,(HPEN )GetStockObject(WHITE_PEN));
	SelectObject(win3.hdc,(HBRUSH )GetStockObject(WHITE_BRUSH));
	Rectangle(win3.hdc,win3.rect.left,win3.rect.top,win3.rect.right,win3.rect.bottom); //Clean drawing area.

	UINT charmode=GetTextAlign(win3.hdc);
	POINT m;

	WCHAR str[100]={0};

	swprintf_s(str,100,L" (%g,%g,%g)",win3.x1,win3.y1,win3.z1);
	SetTextAlign(win3.hdc,TA_RIGHT|TA_TOP);
	moveto3(win3.x1,win3.y1,win3.z1);
	m=getxy3();
	TextOutW(win3.hdc,m.x,m.y,(LPCWSTR)str,wcslen(str));

	moveto3((win3.x1+win3.x2)/2,win3.y1,win3.z1);
	m=getxy3();
	swprintf_s(str,100,L"%S",xstr);
	TextOutW(win3.hdc,m.x,m.y,(LPCWSTR)str,wcslen(str));

	moveto3(win3.x1,win3.y1,(win3.z1+win3.z2)/2);
	m=getxy3();
	swprintf_s(str,100,L"%S",zstr);
	TextOutW(win3.hdc,m.x,m.y,(LPCWSTR)str,wcslen(str));

	SetTextAlign(win3.hdc,TA_LEFT|TA_TOP);
	moveto3(win3.x2,(win3.y1+win3.y2)/2,win3.z1);
	m=getxy3();
	swprintf_s(str,100,L"%S",ystr);
	TextOutW(win3.hdc,m.x,m.y,(LPCWSTR)str,wcslen(str));

	SetTextAlign(win3.hdc,TA_LEFT|TA_BOTTOM);
	swprintf_s(str,100,L"(%g,%g,%g)",win3.x2,win3.y2,win3.z2);
	moveto3(win3.x2,win3.y2,win3.z2);
	m=getxy3();
	TextOutW(win3.hdc,m.x,m.y,(LPCWSTR)str,wcslen(str));

	SetTextAlign(win3.hdc,charmode);
	/////////
	HPEN pen=CreatePen(PS_SOLID,1,dotLineColor);
	HPEN oldpen=(HPEN)SelectObject(win3.hdc,pen);

	line3(win3.x1,win3.y1,win3.z1,win3.x2,win3.y1,win3.z1);
	lineto3(win3.x2,win3.y2,win3.z1);
	lineto3(win3.x2,win3.y2,win3.z2);
	lineto3(win3.x1,win3.y2,win3.z2);
	lineto3(win3.x1,win3.y1,win3.z2);
	lineto3(win3.x1,win3.y1,win3.z1);

	line3(win3.x1,win3.y1,win3.z1,win3.x1,win3.y2,win3.z1);
	lineto3(win3.x2,win3.y2,win3.z1);
	line3(win3.x1,win3.y2,win3.z1,win3.x1,win3.y2,win3.z2);

	double grid,g;
	grid=(win3.x2-win3.x1)/10;
	for(g=0;g<win3.x2;g=g+grid) {
		line3(g,win3.y1,win3.z1,g,win3.y2,win3.z1);
		lineto3(g,win3.y2,win3.z2);
	}
	for(g=0;g>win3.x1;g=g-grid) {
		line3(g,win3.y1,win3.z1,g,win3.y2,win3.z1);
		lineto3(g,win3.y2,win3.z2);
	}

	grid=(win3.y2-win3.y1)/10;
	for(g=0;g<win3.y2;g=g+grid) {
		line3(win3.x2,g,win3.z1,win3.x1,g,win3.z1);
		lineto3(win3.x1,g,win3.z2);
	}
	for(g=0;g>win3.y1;g=g-grid) {
		line3(win3.x2,g,win3.z1,win3.x1,g,win3.z1);
		lineto3(win3.x1,g,win3.z2);
	}

	grid=(win3.z2-win3.z1)/10;
	for(g=0;g<win3.z2;g=g+grid) {
		line3(win3.x1,win3.y1,g,win3.x1,win3.y2,g);
		lineto3(win3.x2,win3.y2,g);
	}
	for(g=0;g>win3.z1;g=g-grid) {
		line3(win3.x1,win3.y1,g,win3.x1,win3.y2,g);
		lineto3(win3.x2,win3.y2,g);
	}

	SelectObject(win3.hdc,oldpen);
	DeleteObject(pen);
	sprintf_s(win3.xstr,50,xstr);
	sprintf_s(win3.ystr,50,ystr);
	sprintf_s(win3.zstr,50,zstr);
}


void frame3(char *xstr,char *ystr,char *zstr) {
	frame3(WHITE,BLACK,xstr,ystr,zstr);
}


void frame3() {
	frame3(win3.xstr,win3.ystr,win3.zstr);
}


void arrow3(COLORREF linecolor,double ax,double ay,double az,double bx,double by,double bz,int h,int d2) {
	POINT a,b;
	moveto3(ax,ay,az);
	a=getxy3();
	moveto3(bx,by,bz);
	b=getxy3();
	arrowline(win3.hdc,3,RGB(255,255,255),linecolor,a,b);
}


void xyz3(COLORREF axiesColor) {
	int h=10,d2=2;

	arrow3(axiesColor,win3.x1,0,0,win3.x2,0,0,h,d2);
	arrow3(axiesColor,0,win3.y1,0,0,win3.y2,0,h,d2);
	arrow3(axiesColor,0,0,win3.z1,0,0,win3.z2,h,d2);
}


void initgraph00(WCHAR title[],int x,int y,int width,int higth) {
	system("color f0"); //Set the console window background to white and foreground to black.

	CONSOLE_CURSOR_INFO cursor_info = {1, 0}; 
	/*Hide the console cursor to avoid disturbing the image*/
	SetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &cursor_info); 
	HWND hWnd=GetConsoleWindow();
	HDC hDC=GetDC(hWnd); //ͨGet the device field handle of the window through the window handle

	SetConsoleTitle((LPCTSTR)title); //Edit window title.
	SetWindowPos(hWnd,HWND_TOPMOST ,x,y,width,higth,SWP_SHOWWINDOW|SWP_NOCOPYBITS); 
	SetWindowLong(hWnd, GWL_STYLE, WS_POPUP); //Edit window style

	RECT r;
	win2.hwnd=hWnd;
	win2.hdc=hDC;
	win3.hwnd=hWnd;
	win3.hdc=hDC;

	GetClientRect(hWnd,&r);
	win2.rect=r;
	win3.rect=r;
}


void initgraph(WCHAR title[],int x,int y,int width,int higth) {
	// system("color f0"); //Set the console window background to white and foreground to black.

	CONSOLE_CURSOR_INFO cursor_info = {1, 0}; 
	/*Hide the console cursor to avoid disturbing the image*/
	SetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &cursor_info);
	HWND hWnd=GetConsoleWindow();
	HDC hDC=GetDC(hWnd); //ͨGet the device field handle of the window through the window handle

	SetConsoleTitle((LPCTSTR)title); //Edit window title.
	SetWindowPos(hWnd,HWND_TOPMOST,x,y,width,higth,SWP_SHOWWINDOW|SWP_NOCOPYBITS); 
	
	RECT r;
	win2.hwnd=hWnd;
	win2.hdc=hDC;
	win3.hwnd=hWnd;
	win3.hdc=hDC;

	GetClientRect(hWnd,&r);
	win2.rect=r;
	win3.rect=r;
}


void plot3(COLORREF color, double (*fxy)(double, double)) {
	double x,y;
	double grid=win3.gridPixels*win3.dyPerPixel;
	HPEN hPen=CreatePen(PS_SOLID, 1, color);
	HPEN hOldPen=(HPEN)SelectObject(win2.hdc, hPen);
	for (y=win3.y1;y<=win3.y2;y=y+grid) {
		x=win3.x1, moveto3(x,y,(*fxy)(x,y));
		for(x=win3.x1;x<=win3.x2;x=x+win3.dxPerPixel) {
			lineto3(x,y,(*fxy)(x,y));
		}
	}
	y=win3.y2;
	{
		x=win3.x1,moveto3(x,y,(*fxy)(x,y));\
		for(x=win3.x1;x<=win3.x2;x=x+win3.dxPerPixel)\
		lineto3(x,y,(*fxy)(x,y));
	}
	grid=win3.gridPixels*win3.dxPerPixel;
	for(x=win3.x1;x<=win3.x2;x=x+grid) {
		y=win3.y1,moveto3(x,y,(*fxy)(x,y));
		for(y=win3.y1;y<=win3.y2;y=y+win3.dyPerPixel) {
			lineto3(x,y,(*fxy)(x,y));
		}	
	}
	x=win3.x2;
	{
		y=win3.y1,moveto3(x,y,(*fxy)(x,y));
		for(y=win3.y1;y<=win3.y2;y=y+win3.dyPerPixel)
		lineto3(x,y,(*fxy)(x,y));
	}
	DeleteObject(hPen);
}


void plot3(COLORREF gridcolor,COLORREF linecolor,COMPLEX p[],int M) {
	int i;
	HPEN hPen=CreatePen(PS_SOLID, 1, gridcolor);
	HPEN hOldPen=(HPEN)SelectObject(win2.hdc, hPen);
	for(i=0;i<M;i++) {
		line3(i,0,0,i,p[i].r,p[i].i);
	}

	HPEN hPen2=CreatePen(PS_SOLID, 1, linecolor);
	SelectObject(win2.hdc, hPen2);
	moveto3(0,p[0].r,p[0].i);
	for(i=0;i<M;i++) {
		lineto3(i,p[i].r,p[i].i);
	}

	SelectObject(win2.hdc,hOldPen);
	DeleteObject(hPen2);
	DeleteObject(hPen);
}


winx::winx(USEWIN win, char *title, int x0, int y0, int width, int height){
	hinstance=(HINSTANCE)GetModuleHandle(NULL);
	if(win==USE_CONSOLE){
		system("color f0"); //Set the console window background to white and foreground to black.

		CONSOLE_CURSOR_INFO cursor_info = {1, 0};
		/*Hide the console cursor to avoid disturbing the image*/
		SetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &cursor_info);

		hwnd0=GetConsoleWindow();
		hdc0=GetDC(hwnd0); //ͨGet the device field handle of the window through the window handle

		SetConsoleTitle((LPCTSTR)title); //Edit window title.
		SetWindowPos(hwnd0,HWND_TOPMOST,x0,y0,width,height,SWP_SHOWWINDOW);
	}
}


void xy2(COLORREF color) {
	arrowline(win2.hdc,3,RGB(255,255,255),color,win2.rect2.left,win2.yo,win2.rect2.right-12,win2.yo);
	arrowline(win2.hdc,3,RGB(255,255,255),color,win2.xo,win2.rect2.bottom,win2.xo,win2.rect2.top+12);
}


void window2(WCHAR title[],double x1,double y1,double x2,double y2,char *xAxis,char *yAxis) {
	initgraph(title,0,0,800,600);
	coord2(x1,y1,x2,y2);
	frame2(xAxis,yAxis);
}


void clearviewport(COLORREF bkcolor) {
	POINT v[4];
	v[0].x=win2.rect2.left,v[0].y=win2.rect2.top;

	v[1].x=win2.rect2.left,v[1].y=win2.rect2.bottom;
	v[2].x=win2.rect2.right,v[2].y=win2.rect2.bottom;
	v[3].x=win2.rect2.right,v[3].y=win2.rect2.top;

	HBRUSH brush=CreateSolidBrush(bkcolor);
	HPEN pen=CreatePen(PS_NULL,1,bkcolor);

	HBRUSH oldbrush=(HBRUSH)SelectObject(win2.hdc,brush);
	HPEN oldpen=(HPEN)SelectObject(win2.hdc,pen);
	Polygon(win2.hdc,v,4);
	SelectObject(win2.hdc,oldpen);DeleteObject(pen);
	SelectObject(win2.hdc,oldbrush);DeleteObject(brush);
}


void setvp3(int scrx,int scry,int width,int heigth,WCHAR title[],char *notice) {
	initgraph(title,scrx,scry,width,heigth);
}


void window3(WCHAR title[],double x1,double y1,double z1,double x2,double y2,double z2,char *xAxis,char *yAxis,char *zAxis) {
	initgraph(title,0,0,800,600);
	coord3(x1,y1,z1,x2,y2,z2);
	frame3(LIGHTGRAY,BLACK,xAxis,yAxis,zAxis);
}


void clear3() {
	SelectObject(win3.hdc,(HPEN )GetStockObject(BLACK_PEN));
	SelectObject(win3.hdc,(HBRUSH )GetStockObject(WHITE_BRUSH));
	Rectangle(win3.hdc,win3.rect.left,win3.rect.top,win3.rect.right,win3.rect.bottom);
}


BOOL ERASEDRAW2=true,ERASEDRAW3=true;


void instKeyCtr3() {
	int N=10;
	int x1=(win3.rect.left+win3.rect.right)/2,y=win3.rect.bottom-N,x2;
	x2=x1+50;
	x1=x1-50;

	HPEN pen1=CreatePen(PS_SOLID,11,BLACK),oldpen=(HPEN)SelectObject(win3.hdc,pen1);
	MoveToEx(win3.hdc,x1,y,NULL);
	LineTo(win3.hdc,x2,y);

	HPEN pen2=CreatePen(PS_SOLID,9,WHITE);SelectObject(win3.hdc,pen2);
	MoveToEx(win3.hdc,x1,y,NULL);
	LineTo(win3.hdc,x2,y);

	HPEN pen3=CreatePen(PS_SOLID,1,BLACK);SelectObject(win3.hdc,pen3);
	MoveToEx(win3.hdc,x1,y,NULL);
	LineTo(win3.hdc,x2,y);

	HPEN pen4=CreatePen(PS_SOLID,1,RED);SelectObject(win3.hdc,pen4);
	MoveToEx(win3.hdc,x1,y-5,NULL);
	LineTo(win3.hdc,x1,y+5);

	SelectObject(win3.hdc,oldpen); //Restore the original pen.
	DeleteObject(pen1);
	DeleteObject(pen2);
	DeleteObject(pen3);
	DeleteObject(pen4);
}


void instKeyCtr2() {
	int N=10;
	int x1=(win2.rect.left+win2.rect.right)/2,y=win2.rect.bottom-N,x2;
	x2=x1+50;
	x1=x1-50;

	HPEN pen1=CreatePen(PS_SOLID,11,BLACK),oldpen=(HPEN)SelectObject(win2.hdc,pen1);
	MoveToEx(win2.hdc,x1,y,NULL);
	LineTo(win2.hdc,x2,y);

	HPEN pen2=CreatePen(PS_SOLID,9,WHITE);SelectObject(win2.hdc,pen2);
	MoveToEx(win2.hdc,x1,y,NULL);
	LineTo(win2.hdc,x2,y);

	HPEN pen3=CreatePen(PS_SOLID,1,BLACK);SelectObject(win3.hdc,pen3);
	MoveToEx(win2.hdc,x1,y,NULL);
	LineTo(win2.hdc,x2,y);

	HPEN pen4=CreatePen(PS_SOLID,1,RED);SelectObject(win3.hdc,pen4);
	MoveToEx(win2.hdc,x1,y-5,NULL);
	LineTo(win2.hdc,x1,y+5);

	SelectObject(win2.hdc,oldpen); //Restore the original pen.
	DeleteObject(pen1);
	DeleteObject(pen2);
	DeleteObject(pen3);
	DeleteObject(pen4);
}


BOOL keyCtr2(double x[1],double a,double b,int n)
{
	int k = 0;
	double dx;
	if(n == 0 || b == a) {
		printf("keyCtr()输入参数错误！\n");
		exit(1);
	}
	dx=(b-a)/n;
	for(; k!=27; k=_getch()) {
		if (k=='+'||k=='='||k=='.'||k=='>') {
			x[0]+=dx;
			x[0]=(x[0]>b)?b:x[0];
		}
		else if (k=='-'||k=='_'||k==','||k=='<') {
			x[0]-=dx;
			x[0]=(x[0]<a)?a:x[0];
		}
		else {
			continue;
		}

		if(ERASEDRAW2) {
			frame2(win2.xstr,win2.ystr);
		}

		int N=10;
		int x1=(win2.rect.left+win2.rect.right)/2,y=win2.rect.bottom-N,x2;
		x2=x1+50;
		x1=x1-50;

		HPEN pen1=CreatePen(PS_SOLID,11,BLACK),oldpen=(HPEN)SelectObject(win2.hdc,pen1);
		MoveToEx(win2.hdc,x1,y,NULL);
		LineTo(win2.hdc,x2,y);

		HPEN pen2=CreatePen(PS_SOLID,9,WHITE);SelectObject(win2.hdc,pen2);
		MoveToEx(win2.hdc,x1,y,NULL);
		LineTo(win2.hdc,x2,y);

		HPEN pen3=CreatePen(PS_SOLID,1,BLACK);SelectObject(win3.hdc,pen3);
		MoveToEx(win2.hdc,x1,y,NULL);
		LineTo(win2.hdc,x2,y);

		HPEN pen4=CreatePen(PS_SOLID,1,RED);SelectObject(win3.hdc,pen4);
		MoveToEx(win2.hdc,x1+(int)((x[0]-a)/dx*100/n+0.5),y-5,NULL);
		LineTo(win2.hdc,x1+(int)((x[0]-a)/dx*100/n+0.5),y+5);

		SelectObject(win2.hdc,oldpen); //Restore the original pen.
		DeleteObject(pen1);
		DeleteObject(pen2);
		DeleteObject(pen3);
		DeleteObject(pen4);
		return true;
	}	
	return false;
}


BOOL keyCtr3(double x[1],double a,double b,int n) {
	int k=0;double dx;
	if(n==0||b==a){printf("keyCtr()输入参数错误！\n");exit(1);}
	dx=(b-a)/n;
	for(;k!=27;k=_getch()) {
		if(k=='+'||k=='='||k=='.'||k=='>') {
			x[0]+=dx;
			x[0]=(x[0]>b)?b:x[0];
		}
		else if(k=='-'||k=='_'||k==','||k=='<') {
			x[0]-=dx;
			x[0]=(x[0]<a)?a:x[0];
		}
		else {
			continue;
		}

		if(ERASEDRAW2) {
			frame3(win3.xstr,win3.ystr,win3.zstr);
		}

		int N=10;
		int x1=(win3.rect.left+win3.rect.right)/2,y=win3.rect.bottom-N,x2;
		x2=x1+50;
		x1=x1-50;

		HPEN pen1=CreatePen(PS_SOLID,11,BLACK),oldpen=(HPEN)SelectObject(win3.hdc,pen1);
		MoveToEx(win3.hdc,x1,y,NULL);
		LineTo(win3.hdc,x2,y);

		HPEN pen2=CreatePen(PS_SOLID,9,WHITE);SelectObject(win3.hdc,pen2);
		MoveToEx(win3.hdc,x1,y,NULL);
		LineTo(win3.hdc,x2,y);

		HPEN pen3=CreatePen(PS_SOLID,1,BLACK);SelectObject(win3.hdc,pen3);
		MoveToEx(win3.hdc,x1,y,NULL);
		LineTo(win3.hdc,x2,y);

		HPEN pen4=CreatePen(PS_SOLID,1,RED);SelectObject(win3.hdc,pen4);
		MoveToEx(win3.hdc,x1+(int)((x[0]-a)/dx*100/n+0.5),y-5,NULL);
		LineTo(win3.hdc,x1+(int)((x[0]-a)/dx*100/n+0.5),y+5);

		SelectObject(win3.hdc,oldpen); //Restore the original pen.
		DeleteObject(pen1);
		DeleteObject(pen2);
		DeleteObject(pen3);
		DeleteObject(pen4);
		return true;
	}	
	return false;
}


void clearvp2() {
	HBRUSH oldbrush=(HBRUSH)SelectObject(win2.hdc,(HBRUSH )GetStockObject(WHITE_BRUSH));
	HPEN oldpen=(HPEN)SelectObject(win2.hdc,(HPEN )GetStockObject(BLACK_PEN));
	Rectangle(win2.hdc,win2.rect2.left-1,win2.rect2.top-1,win2.rect2.right+1,win2.rect2.bottom+1);
	SelectObject(win2.hdc,oldpen);
	SelectObject(win2.hdc,oldbrush);
}


void instKeyCtr(double x[1],double a,double b,int N) {
	int H=10;
	int x1=(win3.rect.left+win3.rect.right)/2,y=win3.rect.bottom-H,x2;
	x2=x1+50;
	x1=x1-50;

	HPEN pen1=CreatePen(PS_SOLID,11,BLACK),oldpen=(HPEN)SelectObject(win3.hdc,pen1);
	MoveToEx(win3.hdc,x1,y,NULL);
	LineTo(win3.hdc,x2,y);

	HPEN pen2=CreatePen(PS_SOLID,9,WHITE);SelectObject(win3.hdc,pen2);
	MoveToEx(win3.hdc,x1,y,NULL);
	LineTo(win3.hdc,x2,y);

	HPEN pen3=CreatePen(PS_SOLID,1,BLACK);SelectObject(win3.hdc,pen3);
	MoveToEx(win3.hdc,x1,y,NULL);
	LineTo(win3.hdc,x2,y);

	int gg=(int)((*x-a)/(b-a)*100);
	HPEN pen4=CreatePen(PS_SOLID,1,RED);SelectObject(win3.hdc,pen4);
	MoveToEx(win3.hdc,x1+gg,y-5,NULL);
	LineTo(win3.hdc,x1+gg,y+5); //Draw a red position marker line

	SelectObject(win3.hdc,oldpen); //Restore the original pen.
	DeleteObject(pen1);
	DeleteObject(pen2);
	DeleteObject(pen3);
	DeleteObject(pen4);
}


BOOL keyCtr(double x[1],double a,double b,int n) {
	int k=0;double dx;
	if(n==0||b==a) {
		printf("keyCtr()输入参数错误！\n");
		exit(1);
	}
	dx=(b-a)/n;
	for(;k!=27;k=_getch()) {
		if(k=='+'||k=='='||k=='.'||k=='>') {
			x[0]+=dx;
			x[0]=(x[0]>b)?b:x[0];
		}
		else if(k=='-'||k=='_'||k==','||k=='<') {
			x[0]-=dx;
			x[0]=(x[0]<a)?a:x[0];
		}
		else {
			continue;
		} 
		return true;
	}	
	return false;
}


void plotgri2(COLORREF gridcolor, COLORREF linecolor, double *p, int N) {
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
}


void plot_frequency_spectrum(COLORREF gridcolor, double *p, int N) {
	HPEN pen = CreatePen(PS_SOLID,1,gridcolor), oldpen=(HPEN)SelectObject(win3.hdc,pen);
	for(int i = 0; i < N; i++) {
		line2(i, 0.0, i, p[i]);
	}
}