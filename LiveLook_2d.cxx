// This is a file to read the data outize
//ut from a the TDC

// == INCLUDES ==
#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <dirent.h>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath> 
#include <iostream>
#include <ctype.h>

#define TDCsignalLen 15
#define READ_LIMIT 10000
#define EVENT_LIMIT 100000
#define cursup "\033[A"

class HistoGUI{


	public:

	HistoGUI(){}

	Display * disp;
	Window    wind;
	XEvent    evt;
	int       screen;

	XColor xcolour;
	XColor xcolour_red;
	XColor xcolour_veryred;
	Colormap cmap;

	XColor PixelColour[10];

	std::vector<double> x;
	std::vector<double> y;

	std::vector<double> x_2;
	std::vector<double> y_2;
	std::vector<std::vector<double>> content_2;

	std::vector<double> x_3;
	std::vector<double> y_3;
	std::vector<std::vector<double>> content_3;

	double max_x;
	double min_x;
	double max_y;
	double min_y;

	double pos_x;
	double pos_y;
	double new_pos_x;
	double new_pos_y;

	unsigned int width;  
	unsigned int height;  
	double window_width;
    double window_height;
	double width_scale; 
	double x_offset;
	double height_scale;
	double y_offset;
	bool Draw2D_On;

	double old_xl, old_xh, old_yl, old_yh;
	double old_mouse_x, old_mouse_y;

	int Init();
	int SetData(std::vector<double> a, std::vector<double> b);
	int SetData2(std::vector<double> a, std::vector<double> b, std::vector<std::vector<double>> c);
	int Loop();	
	void Close(){ printf("Closing now!\n"); XCloseDisplay(disp); }
	int DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
	int DrawData2D(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
	int SwapData();

	int DrawCrosshairs(int mouse_x, int mouse_y);
	int Zoom(int mouse_x, int mouse_y);
	long int ReReadChannel(char* filename, char* folderName, long int seek_to);
	int TDCLoop(long int *locs, char filenames[16][64], char * folderName);

    std::map<int, int> pixelmap = {
        //ChNo  PIXEL
		// H - Type
		{ 0 , 52 ,},
		{ 1 , 57 ,},
		{ 2 , 60 ,},
		{ 3 , 58 ,},
		{ 4 , 59 ,},
		{ 5 , 50 ,}, 
		{ 6 , 49 ,}, 
		{ 7 , 51 ,},
		{ 8 , 36 ,},  
		{ 9 , 41 ,},
		{ 10, 44 ,}, 
		{ 11, 42 ,},
		{ 12, 43 ,},
		{ 13, 33 ,},
		{ 14, 27 ,},
		{ 15, 34 ,},
		{ 16, 35 ,},
		{ 17, 18 ,},
		{ 18, 19 ,},
		{ 19, 25 ,},
		{ 20, 26 ,},
		{ 21, 28 ,},
		{ 22, 9  ,},
		{ 23, 10 ,},
		{ 24, 20 ,},
		{ 25, 1  ,},
		{ 26, 3  ,},
		{ 27, 12 ,},
		{ 28, 11 ,},
		{ 29, 4  ,},
		{ 30, 2  ,},
		{ 31, 17 ,}
	};


};

int HistoGUI::SetData(std::vector<double> a, std::vector<double>b){
	for(int i=0; i < a.size(); i++){
		x.push_back(a[i]);
		y.push_back(b[i]);
	}
	return a.size();
}

int HistoGUI::SetData2(std::vector<double> a, std::vector<double> b, std::vector<std::vector<double>> c){
	for(int i=0; i < b.size(); i++){
		y_2.push_back(b[i]);
	}
	for(int j=0; j < a.size(); j++){
		x_2.push_back(a[j]);
		std::vector<double> temp;
		for(int i=0; i < b.size(); i++){
			temp.push_back(c[j][i]);
		}
		content_2.push_back(temp);
	}
	return a.size();
}

int HistoGUI::SwapData(){
	std::vector<double> temp_x = x_2;	
	std::vector<double> temp_y = y_2;	
	std::vector<std::vector<double>> temp_cont = content_2;	

	content_2 = content_3;
	x_2 = x_3;
	y_2 = y_3;
	content_3 = temp_cont;
	x_3 = temp_x;
	y_3 = temp_y;

	return x_2.size();
}


int HistoGUI::Init(){

	Draw2D_On = false;
	disp = XOpenDisplay(NULL);
	if (disp == NULL) {
	   fprintf(stderr, "Cannot open display\n");
	   exit(1);
	}
	
	screen = DefaultScreen(disp);
	wind   = XCreateSimpleWindow(disp, RootWindow(disp, screen), 10, 10, 600, 400, 1,
	                        BlackPixel(disp, screen), WhitePixel(disp, screen));
	XGrabPointer(disp, wind, False, ButtonPressMask | ButtonReleaseMask | PointerMotionMask, GrabModeAsync,
                 GrabModeAsync, None, None, CurrentTime);
	XSelectInput(disp, wind, ExposureMask | KeyPressMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask);
	XMapWindow(disp, wind);

	// colours
    cmap = DefaultColormap(disp, screen);
	xcolour.red = 32000; 
	xcolour.green = 32000; 
	xcolour.blue = 42000;
	xcolour.flags = DoRed | DoGreen | DoBlue;
	XAllocColor(disp, cmap, &xcolour);

	xcolour_red.red   = 42000; 
	xcolour_red.green = 32000; 
	xcolour_red.blue  = 32000; 
	xcolour_red.flags = DoRed | DoGreen | DoBlue;
	XAllocColor(disp, cmap, &xcolour_red);

	xcolour_veryred.red   = 65000; 
	xcolour_veryred.green = 32000; 
	xcolour_veryred.blue  = 32000; 
	xcolour_veryred.flags = DoRed | DoGreen | DoBlue;
	XAllocColor(disp, cmap, &xcolour_veryred);

	for(int i = 0; i < 10; i++){	
		PixelColour[i].red   = 65535; 
		PixelColour[i].green = 65535 * (10 - i)/ 10; 
		PixelColour[i].blue  = 65535 * (10 - i)/ 10; 
		PixelColour[i].flags = DoRed | DoGreen | DoBlue;
		XAllocColor(disp, cmap, &PixelColour[i]);
	}


	return 1;

}

int HistoGUI::DrawCrosshairs(int mouse_x, int mouse_y){

	int j1, j2;
	unsigned int j3, j4;  
	Window root_return;

	XGetGeometry(disp, wind, &root_return, &j1, &j2, &width, &height, &j3, &j4);
	XSetForeground(disp,  DefaultGC(disp, screen), xcolour_red.pixel);
	XDrawLine(disp, wind, DefaultGC(disp, screen), mouse_x, 0.0, mouse_x, height);
	XDrawLine(disp, wind, DefaultGC(disp, screen), 0.0, mouse_y, width,  mouse_y);
	
	pos_y = mouse_y * height_scale - y_offset;
	pos_x = mouse_x * width_scale  - x_offset;

	char coord[32];
	sprintf(coord, "(%.2f, %.2f)", pos_x, pos_y);
	XDrawString(disp, wind, DefaultGC(disp, screen), mouse_x + 5, mouse_y + 12, coord, strlen(coord));
	 
	return 1;
}

int HistoGUI::Zoom(int mouse_x, int mouse_y){
	
	new_pos_x = mouse_x * width_scale - x_offset;
	new_pos_y = mouse_y * height_scale  - y_offset;

	pos_x = old_mouse_x * width_scale - x_offset;
	pos_y = old_mouse_y * height_scale  - y_offset;
	

	if(std::abs(mouse_x - old_mouse_x) < 20 or std::abs(mouse_y - old_mouse_y) < 20) return 1;


	//printf("Zooming: %f - %f, %f - %f\n", mouse_x, old_mouse_x,mouse_y,old_mouse_y);

	double x_low_win, y_low_win, x_hi_win, y_hi_win;
	if(new_pos_x < pos_x){
		x_low_win = new_pos_x;
		x_hi_win  = pos_x;
		//printf("new low; [%f, %f]\n",x_low_win, x_hi_win);
	} else {
		x_low_win = pos_x;
		x_hi_win  = new_pos_x;
		//printf("new high; [%f, %f]\n",x_low_win, x_hi_win);
	}
	if(new_pos_y > pos_y){
		y_low_win = new_pos_y;
		y_hi_win  = pos_y;
	} else {
		y_low_win = pos_y;
		y_hi_win  = new_pos_y;
	}
	

	//printf("[(%f,%f) (%f,%f)]\n",x_low_win, y_low_win, x_hi_win, y_hi_win);
	XClearWindow(disp, wind);
	DrawData(x_low_win, y_low_win, x_hi_win, y_hi_win);
	//double x_low_win, double y_low_win, double x_hi_win, double y_hi_win

	return 1;

}



int HistoGUI::DrawData2D(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win){
	int j1, j2;
	unsigned int j3, j4;  
	Window root_return;
	XGetGeometry(disp, wind, &root_return, &j1, &j2, &width, &height, &j3, &j4);

	double x_step;// = (width  * 0.8) / x.size();
	double y_step;// = (height * 0.8) / x.size();

	if(x_low_win == -1 and y_low_win == -1 and x_hi_win == -1 and y_hi_win == -1){	
		// Audomatically decide data postition
		max_x = x_2[0];
		max_y = y_2[0];
		min_x = x_2[0];
		min_y = y_2[0];

		for(int i=0; i<x_2.size(); i++){
			if(x_2[i] > max_x) max_x = x_2[i];
			if(x_2[i] < min_x) min_x = x_2[i];
		}
		for(int i=0; i<y_2.size(); i++){
			if(y_2[i] > max_y) max_y = y_2[i];
			if(y_2[i] < min_y) min_y = y_2[i];
		}

		double max_cont = content_2[0][0];	

		for(int i=0; i<x_2.size(); i++){
			for(int j=0; j<y_2.size(); j++){
				if(content_2[i][j] > max_cont) max_cont = content_2[i][j];
			}
		}

		//printf(" max_x = %f\n", max_x );
		//printf(" max_y = %f\n", max_y );
		//printf(" min_x = %f\n", min_x );
		//printf(" min_y = %f\n", min_y );

		width_scale = (max_x - min_x) / (0.8 * width);
		x_offset = 0.5 * ((max_x - min_x) / 0.8) - 0.5 * (min_x + max_x);

		height_scale = -1. * (max_y - min_y) / (0.8 * height);
		y_offset = -0.5 * ((max_y - min_y) / 0.8) + -0.5 * (min_y + max_y);
		
		//printf("width scale  = %f\n", width_scale);
		//printf("x_offset     = %f\n", x_offset);
		//printf("height scale = %f\n", height_scale);
		//printf("y_offset     = %f\n", y_offset);	
	
		XSetForeground(disp, DefaultGC(disp,screen), 0);	
		double x_wid ;
		double y_wid ;
		double x_wid2;
		double y_wid2;

		double binwidth_x = 0.8 * width_scale  / x_2.size();
		double binwidth_y = 0.8 * height_scale / y_2.size();
	

		for(int i=0; i < width; i++){
			x_wid  = 1+ (i * width_scale) - x_offset;
			if(x_wid > max_x) x_wid = max_x;
			if(x_wid < min_x) x_wid = min_x;

			for(int j=0; j < height; j++){
				y_wid  = 1 + (j * height_scale) - y_offset;
				if(y_wid > max_y) y_wid = max_y;
				if(y_wid < min_y) y_wid = min_y;
				int colindex = (int) 10 * content_2[(int)x_wid][(int)y_wid] / max_cont; 
				//int colindex = (int) 10 * i / width; 
				//if(colindex > 0){
				//printf("Colour = %i\n",colindex);
				//printf("(%f, %f) = %f\n", x_wid,y_wid, content_2[(int)x_wid][(int)y_wid]);
				//}
				XSetForeground(disp, DefaultGC(disp,screen), PixelColour[colindex].pixel);
				//XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid - 0.5* binwidth_x, y_wid -0.5*binwidth_y, binwidth_x, binwidth_y);
				XDrawPoint(disp, wind, DefaultGC(disp, screen), i,j);
			}
		}

		double axis_x  = (0. + x_offset) / width_scale;
		double axis_y  = (0. + y_offset) / height_scale;

		XSetForeground(disp, DefaultGC(disp,screen), xcolour.pixel);
		XDrawLine(disp, wind, DefaultGC(disp, screen), axis_x, 0.0, axis_x, height);
		XDrawLine(disp, wind, DefaultGC(disp, screen), 0.0, axis_y, width, axis_y);

		char axis_val[4];
		int w_step = width / 10;
		for(int i=0; i < (int) width; i += w_step){
			double x_val   = i * width_scale - x_offset;
			sprintf(axis_val, "%.1f", x_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), i, axis_y + 10, axis_val, strlen(axis_val));
		}

		int h_step = height / 10;
		for(int i=0; i < (int) height; i += h_step){
			double y_val = i * height_scale - y_offset;
			sprintf(axis_val, "%.1f", y_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), axis_x + 10, i, axis_val, strlen(axis_val));
		}

	} else {
		//double x_low_win, double y_low_win, double x_hi_win, double y_hi_win

		double max_cont = content_2[0][0];	
		for(int i=0; i<x_2.size(); i++){
			for(int j=0; j<y_2.size(); j++){
				if(content_2[i][j] > max_cont) max_cont = content_2[i][j];
			}
		}

		width_scale = (x_hi_win - x_low_win) / width;
		x_offset = -1.0 * x_low_win;

		height_scale = (y_hi_win - y_low_win) / height;
		y_offset = -1.0 * y_low_win;
		
		//printf("width scale  = %f\n", width_scale);
		//printf("x_offset     = %f\n", x_offset);
		//printf("height scale = %f\n", height_scale);
		//printf("y_offset     = %f\n", y_offset);	
	
		XSetForeground(disp, DefaultGC(disp,screen), 0);	
		double x_wid ;
		double y_wid ;
		double x_wid2;
		double y_wid2;

		double binwidth_x = 0.8 * width_scale  / x_2.size();
		double binwidth_y = 0.8 * height_scale / y_2.size();

		for(int i=0; i < width; i++){
			x_wid  = (i * width_scale) - x_offset;
			if(x_wid > max_x) x_wid = max_x;
			if(x_wid < min_x) x_wid = min_x;

			for(int j=0; j < height; j++){
				y_wid  = (j * height_scale) - y_offset;
				if(y_wid > max_y) y_wid = max_y;
				if(y_wid < min_y) y_wid = min_y;
				int colindex = (int) 10 * content_2[(int)x_wid][(int)y_wid] / max_cont; 
				//int colindex = (int) 10 * i / width; 
				//if(colindex > 0){
				//printf("Colour = %i\n",colindex);
				//printf("(%f, %f) = %f\n", x_wid,y_wid, content_2[(int)x_wid][(int)y_wid]);
				//}
				XSetForeground(disp, DefaultGC(disp,screen), PixelColour[colindex].pixel);
				//XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid - 0.5* binwidth_x, y_wid -0.5*binwidth_y, binwidth_x, binwidth_y);
				XDrawPoint(disp, wind, DefaultGC(disp, screen), i,j);
			}
		}

		double axis_x  = (0. + x_offset) / width_scale;
		double axis_y  = (0. + y_offset) / height_scale;

		XSetForeground(disp, DefaultGC(disp,screen), xcolour.pixel);
		XDrawLine(disp, wind, DefaultGC(disp, screen), axis_x, 0.0, axis_x, height);
		XDrawLine(disp, wind, DefaultGC(disp, screen), 0.0, axis_y, width, axis_y);

		char axis_val[4];
		int w_step = width / 10;
		for(int i=0; i < (int) width; i += w_step){
			double x_val   = i * width_scale - x_offset;
			sprintf(axis_val, "%.1f", x_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), i, axis_y + 10, axis_val, strlen(axis_val));
		}

		int h_step = height / 10;
		for(int i=0; i < (int) height; i += h_step){
			double y_val = i * height_scale - y_offset;
			sprintf(axis_val, "%.1f", y_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), axis_x + 10, i, axis_val, strlen(axis_val));
		}


	}

	old_xl = x_low_win;
	old_yl = y_low_win;
	old_xh = x_hi_win;
	old_yh = y_hi_win;

	return 1;
}


int HistoGUI::DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win){

	if(Draw2D_On){
		DrawData2D(x_low_win, y_low_win, x_hi_win, y_hi_win);
		return 1;
	}

	int j1, j2;
	unsigned int j3, j4;  
	Window root_return;

	XGetGeometry(disp, wind, &root_return, &j1, &j2, &width, &height, &j3, &j4);

	//printf("Width = %u, Height = %u, x = %i, y = %i\n", width, height, j1, j2);
	//printf("[(%f,%f) (%f,%f)]\n",x_low_win, y_low_win, x_hi_win, y_hi_win);

	double x_step;// = (width  * 0.8) / x.size();
	double y_step;// = (height * 0.8) / x.size();


	if(x_low_win == -1 and y_low_win == -1 and x_hi_win == -1 and y_hi_win == -1){	
		// Audomatically decide data postition
		max_x = x[0];
		max_y = y[0];
		min_x = x[0];
		min_y = y[0];

		for(int i=0; i<x.size(); i++){
			if(x[i] > max_x) max_x = x[i];
			if(x[i] < min_x) min_x = x[i];
			if(y[i] > max_y) max_y = y[i];
			if(y[i] < min_y) min_y = y[i];
		}
		//printf(" max_x = %f\n", max_x );
		//printf(" max_y = %f\n", max_y );
		//printf(" min_x = %f\n", min_x );
		//printf(" min_y = %f\n", min_y );

		width_scale = (max_x - min_x) / (0.8 * width);
		x_offset = 0.5 * ((max_x - min_x) / 0.8) - 0.5 * (min_x + max_x);

		height_scale = -1. * (max_y - min_y) / (0.8 * height);
		y_offset = -0.5 * ((max_y - min_y) / 0.8) + -0.5 * (min_y + max_y);
		
		//printf("width scale  = %f\n", width_scale);
		//printf("x_offset     = %f\n", x_offset);
		//printf("height scale = %f\n", height_scale);
		//printf("y_offset     = %f\n", y_offset);
		double axis_x  = (0. + x_offset) / width_scale;
		double axis_y  = (0. + y_offset) / height_scale;

		XSetForeground(disp, DefaultGC(disp,screen), xcolour.pixel);
		XDrawLine(disp, wind, DefaultGC(disp, screen), axis_x, 0.0, axis_x, height);
		XDrawLine(disp, wind, DefaultGC(disp, screen), 0.0, axis_y, width, axis_y);

		char axis_val[4];
		int w_step = width / 10;
		for(int i=0; i < (int) width; i += w_step){
			double x_val   = i * width_scale - x_offset;
			sprintf(axis_val, "%.1f", x_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), i, axis_y + 10, axis_val, strlen(axis_val));
		}

		int h_step = height / 10;
		for(int i=0; i < (int) height; i += h_step){
			double y_val = i * height_scale - y_offset;
			sprintf(axis_val, "%.1f", y_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), axis_x + 10, i, axis_val, strlen(axis_val));
		}
	
	
		XSetForeground(disp, DefaultGC(disp,screen), 0);	
		double x_wid ;
		double y_wid ;
		double x_wid2;
		double y_wid2;

		for(int i=0; i < x.size() - 2; i++){
			x_wid  = (x[i] + x_offset) / width_scale;
			y_wid  = (y[i] + y_offset) / height_scale;
			x_wid2 = (x[i + 1] + x_offset) / width_scale;
			y_wid2 = (y[i + 1] + y_offset) / height_scale;
			//printf("(%f, %f), (%f,%f)\n", x_wid,y_wid,x_wid2,y_wid2);
			XDrawLine(disp, wind, DefaultGC(disp, screen), x_wid, y_wid, x_wid2, y_wid2);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -2, y_wid -2, 4, 4);
		}
		XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid2 -2, y_wid2 -2, 4, 4);

	} else {
		//double x_low_win, double y_low_win, double x_hi_win, double y_hi_win

		width_scale = (x_hi_win - x_low_win) / width;
		x_offset = -1.0 * x_low_win;

		height_scale = (y_hi_win - y_low_win) / height;
		y_offset = -1.0 * y_low_win;
		
		//printf("width scale  = %f\n", width_scale);
		//printf("x_offset     = %f\n", x_offset);
		//printf("height scale = %f\n", height_scale);
		//printf("y_offset     = %f\n", y_offset);
		double axis_x  = (0. + x_offset) / width_scale;
		double axis_y  = (0. + y_offset) / height_scale;

		XSetForeground(disp, DefaultGC(disp,screen), xcolour.pixel);
		XDrawLine(disp, wind, DefaultGC(disp, screen), axis_x, 0.0, axis_x, height);
		XDrawLine(disp, wind, DefaultGC(disp, screen), 0.0, axis_y, width, axis_y);

		char axis_val[4];
		int w_step = width / 10;
		for(int i=0; i < (int) width; i += w_step){
			double x_val   = i * width_scale - x_offset;
			sprintf(axis_val, "%.1f", x_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), i, axis_y + 10, axis_val, strlen(axis_val));
		}

		int h_step = height / 10;
		for(int i=0; i < (int) height; i += h_step){
			double y_val = i * height_scale - y_offset;
			sprintf(axis_val, "%.1f", y_val);
			XDrawString(disp, wind, DefaultGC(disp, screen), axis_x + 10, i, axis_val, strlen(axis_val));
		}
	
	
		XSetForeground(disp, DefaultGC(disp,screen), 0);	
		double x_wid ;
		double y_wid ;
		double x_wid2;
		double y_wid2;
	
		for(int i=0; i < x.size() - 2; i++){
			x_wid  = (x[i] + x_offset) / width_scale;
			y_wid  = (y[i] + y_offset) / height_scale;
			x_wid2 = (x[i + 1] + x_offset) / width_scale;
			y_wid2 = (y[i + 1] + y_offset) / height_scale;
		//	printf("(%f, %f), (%f,%f)\n", x_wid,y_wid,x_wid2,y_wid2);
			XDrawLine(disp, wind, DefaultGC(disp, screen), x_wid, y_wid, x_wid2, y_wid2);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -2, y_wid -2, 4, 4);
		}
		XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid2 -2, y_wid2 -2, 4, 4);


	}

	old_xl = x_low_win;
	old_yl = y_low_win;
	old_xh = x_hi_win;
	old_yh = y_hi_win;

	return 1;
}


int HistoGUI::Loop(){

	bool MousePressed  = false;
	bool MousePressed2 = false;

	while (1) {
		XNextEvent(disp, &evt);
		if (evt.type == Expose) {
	//		XFillRectangle(disp, wind, DefaultGC(disp, screen), 20, 20, 10, 10);
//			XDrawString   (disp, wind, DefaultGC(disp, screen), 10, 50, msg, strlen(msg));
			DrawData(-1,-1,-1,-1);
		} else if (evt.type == ButtonPress){
			/* store the mouse button coordinates in 'int' variables. */
			/* also store the ID of the window on which the mouse was */
			/* pressed.                                               */
			int mouse_x = evt.xbutton.x;
			int mouse_y = evt.xbutton.y;
			old_mouse_x = mouse_x;
			old_mouse_y = mouse_y;

			/* check which mouse button was pressed, and act accordingly. */
			if(evt.xbutton.button == Button1){
				/* draw a pixel at the mouse position. */
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
				DrawCrosshairs(mouse_x, mouse_y);
				MousePressed = true;
				printf("(%f, %f)\n", pos_x, pos_y);
			} if(evt.xbutton.button == Button3){
				MousePressed2 = true;
			} if(evt.xbutton.button == Button2){
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
			}
		} else if (evt.type == ButtonRelease){
			/* store the mouse button coordinates in 'int' variables. */
			/* also store the ID of the window on which the mouse was */
			/* pressed.                                               */
			int mouse_x = evt.xbutton.x;
			int mouse_y = evt.xbutton.y;
			//printf("Released button: %i, %i\n", mouse_x, mouse_y);
			//printf("Pressed button: %i, %i\n", mouse_x, mouse_y);
			/* check which mouse button was pressed, and act accordingly. */
			if(evt.xbutton.button == Button1){
				/* draw a pixel at the mouse position. */
				Zoom(mouse_x, mouse_y);
				//DrawData();
				MousePressed = false;
			} if(evt.xbutton.button == Button3){
				MousePressed2 = false;
			}
		} else if (evt.type == MotionNotify and MousePressed){

			int mouse_x = evt.xmotion.x;
			int mouse_y = evt.xmotion.y;
			XClearWindow(disp, wind);

			DrawData(old_xl, old_yl, old_xh, old_yh);
			DrawCrosshairs(old_mouse_x, old_mouse_y);
			DrawCrosshairs(mouse_x, mouse_y);

		} else if (evt.type == MotionNotify and MousePressed2){

			int mouse_x = evt.xmotion.x;
			int mouse_y = evt.xmotion.y;

			XSetForeground(disp, DefaultGC(disp,screen), xcolour_veryred.pixel);
		XFillRectangle(disp, wind, DefaultGC(disp, screen), mouse_x -2, mouse_y -2, 4, 4);
	
		} else if (evt.type == KeyPress){

			if(evt.xkey.keycode == 0x41){
				XClearWindow(disp, wind);
				DrawData(-1,-1,-1,-1);
			} else if(evt.xkey.keycode == 0x27){
				SwapData();
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
			} else if(evt.xkey.keycode == 0x19 and !Draw2D_On){
				XClearWindow(disp, wind);
				Draw2D_On = true;
				DrawData2D(-1,-1,-1,-1);
			} else if(evt.xkey.keycode == 0x19 and Draw2D_On){
				XClearWindow(disp, wind);
				Draw2D_On = false;
				DrawData(-1,-1,-1,-1);	
			} else {
				break;	
			}
		}
	}

	return 1;
}

int HistoGUI::TDCLoop(long int *locs, char filenames[16][64], char * folderName){
	printf("Started loop\n");
	fflush(stdout);

	bool MousePressed  = false;
	bool MousePressed2 = false;

	while (1) {
		XNextEvent(disp, &evt);
		if (evt.type == Expose) {
	//		XFillRectangle(disp, wind, DefaultGC(disp, screen), 20, 20, 10, 10);
//			XDrawString   (disp, wind, DefaultGC(disp, screen), 10, 50, msg, strlen(msg));
			DrawData(-1,-1,-1,-1);
		} else if (evt.type == ButtonPress){
			/* store the mouse button coordinates in 'int' variables. */
			/* also store the ID of the window on which the mouse was */
			/* pressed.                                               */
			int mouse_x = evt.xbutton.x;
			int mouse_y = evt.xbutton.y;
			old_mouse_x = mouse_x;
			old_mouse_y = mouse_y;

			/* check which mouse button was pressed, and act accordingly. */
			if(evt.xbutton.button == Button1){
				/* draw a pixel at the mouse position. */
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
				DrawCrosshairs(mouse_x, mouse_y);
				MousePressed = true;
				printf("(%f, %f)\n", pos_x, pos_y);
			} if(evt.xbutton.button == Button3){
				MousePressed2 = true;
			} if(evt.xbutton.button == Button2){
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
			}
		} else if (evt.type == ButtonRelease){
			/* store the mouse button coordinates in 'int' variables. */
			/* also store the ID of the window on which the mouse was */
			/* pressed.                                               */
			int mouse_x = evt.xbutton.x;
			int mouse_y = evt.xbutton.y;
			//printf("Released button: %i, %i\n", mouse_x, mouse_y);
			//printf("Pressed button: %i, %i\n", mouse_x, mouse_y);
			/* check which mouse button was pressed, and act accordingly. */
			if(evt.xbutton.button == Button1){
				/* draw a pixel at the mouse position. */
				Zoom(mouse_x, mouse_y);
				//DrawData();
				MousePressed = false;
			} if(evt.xbutton.button == Button3){
				MousePressed2 = false;
			}
		} else if (evt.type == MotionNotify and MousePressed){

			int mouse_x = evt.xmotion.x;
			int mouse_y = evt.xmotion.y;
			XClearWindow(disp, wind);

			DrawData(old_xl, old_yl, old_xh, old_yh);
			DrawCrosshairs(old_mouse_x, old_mouse_y);
			DrawCrosshairs(mouse_x, mouse_y);

		} else if (evt.type == MotionNotify and MousePressed2){

			int mouse_x = evt.xmotion.x;
			int mouse_y = evt.xmotion.y;

			XSetForeground(disp, DefaultGC(disp,screen), xcolour_veryred.pixel);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), mouse_x -2, mouse_y -2, 4, 4);
	
		} else if (evt.type == KeyPress){
			printf("Key = %x\n", evt.xkey.keycode);

			if(evt.xkey.keycode == 0x41){
				XClearWindow(disp, wind);
				DrawData(-1,-1,-1,-1);
			} else if(evt.xkey.keycode == 0x27){
				SwapData();
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
			} else if(evt.xkey.keycode == 0x19 and !Draw2D_On){
				XClearWindow(disp, wind);
				Draw2D_On = true;
				DrawData2D(-1,-1,-1,-1);
			} else if(evt.xkey.keycode == 0x19 and Draw2D_On){
				XClearWindow(disp, wind);
				Draw2D_On = false;
				DrawData(-1,-1,-1,-1);	
				
			} else if(evt.xkey.keycode == 0x18) {
				break;	
			} else {
				printf("Going\n");
				fflush(stdout);
				for(int file=0; file < 16; file++){
					printf("%s : %f = %li\n",filenames[file],folderName, locs[file]);
					locs[file] = ReReadChannel(filenames[file], folderName, locs[file]); 					
				}
				XClearWindow(disp, wind);
				DrawData(old_xl, old_yl, old_xh, old_yh);
			}
		}
	}

	return 1;
}


int number_of_zero_len_payloads = 0;

// ======== ====== ========
// ======== HEADER ========
// ======== ====== ========

class HeaderTDC{

	public:

	// Public variables
	int ChNo;

	
	// Constructor & destructor events
	HeaderTDC(int a){
			ChNo = a;
	}
//	~HeaderTDC();

	// Read data from a file
	int Read(FILE *infile){
		
		int BR = fread(&header_data, sizeof(unsigned char), 4, infile);
		return BR;
	}

	// Print data to screen
	int Print(){	
		printf("Ch %i header: %.2x%.2x %.2x%.2x\n", ChNo, header_data[0], header_data[1], header_data[2], header_data[3] );
		return 1;
	}

		
	private:
	// Private variables
	unsigned char header_data[4];

};




// ======== ===== ========
// ======== EVENT ========
// ======== ===== ========

// This is the event class, it holds an event.
class EventTDC{

	public:

	// Constructor & destructor events
	EventTDC(int a){
		id = a;
		verbose = 0;
		payloadLength = 0;
		noTriggers = 0;
	}
	EventTDC(int a, int b){
		id = a;
		channel_id = b;
		verbose = 0;
		Ch1_ready = 1;
		Ch2_ready = 1;
		payloadLength = 0;
		noTriggers = 0;
	}
	EventTDC(){
		id = 1111;
		channel_id = 1111;
		verbose = 0;
		Ch1_ready = 1;
		Ch2_ready = 1;
		payloadLength = 0;
		noTriggers = 0;
	}
//	~EventTDC();

	// Public functions
	int Read(FILE *infile);
	int Print();
	int SetVerbose(int ch);
	int size(int ch);
	int GetID(){ return id;};
	int GetChID(){ return channel_id;};
	int SetChID(int a){channel_id = a; return channel_id;};
	int IsVerbose(){ return verbose;};
	int SetChIDFromFilename(char* filename);

	// Data vectors
	std::vector<int> fineTime_1;
	std::vector<int> fineTime_2;
	std::vector<int> tot_1;
	std::vector<int> tot_2;
	std::vector<unsigned int> trigNo_1;
	std::vector<unsigned int> trigNo_2;
	
	int DataFrame1[24 * TDCsignalLen] = { 0 };
	int DataFrame2[24 * TDCsignalLen] = { 0 };
	int payloadLength = 0;
	int noTriggers;
	int tot_1_candidate = 0;
	int tot_2_candidate = 0;
	unsigned long int id = 0;


	private:

	// Private data 
	int searchTerm_1;
	int delay_1; 
	int searchTerm_2;
	int delay_2; 
	int channel_id;
	int verbose;
	int Ch1_ready;
	int Ch2_ready;

	unsigned char bin_data[16];
	unsigned char channel_data[2];
	unsigned char channel_data_A[2];
	
	// Private functions 
	int GetTDCTest(unsigned char data[4]);
	int GetTDCTime1();
	int GetTDCTime2();


};


int EventTDC::Print(){ 

	printf("\n");
	printf("Event %.4lu Ch no : %.2i - %.2i\n", id, channel_id, channel_id+1); 
	printf("==========================\n");
	printf("Payload length = %i\n", payloadLength); 
	printf("Bytes read     = %i\n", 16 * payloadLength); 
	printf("\n");
	
	printf("Channel %i: %i\n",channel_id, size(1));
	printf("==================\n");
	printf("Evt# | Time | ToT\n");

	// Print channel 1 data - not too much mind
	if(size(1) > 5){
		for(int i=0; i < 5; i++) printf("%.4i | %.4i | %.4i\n", i, fineTime_1[i], tot_1[i]);
		printf(" .  .  .  .  .  .\n");
		printf("%.4i | %.4i | %.4i\n", (size(1) - 1), fineTime_1[size(1) - 1], tot_1[size(1) - 1]);
	} else {
		for(int i=0; i < size(1); i++) printf("%.4i | %.4i | %.4i\n", i, fineTime_1[i], tot_1[i]);
	}
	printf("==================\n");

	printf("\n");
	printf("Channel %i: %i\n", channel_id + 1, size(2));
	printf("==================\n");
	printf("Evt# | Time | ToT\n");

	// Print channel 2 data
	if(size(2) > 5){
		for(int i=0; i < 5; i++) printf("%.4i | %.4i | %.4i\n", i, fineTime_2[i], tot_2[i]);
		printf(" .  .  .  .  .  .\n");
		printf("%.4i | %.4i | %.4i\n", (size(2) - 1), fineTime_2[size(2) - 1], tot_2[size(2) - 1]);
	} else {
		for(int i=0; i < size(2); i++) printf("%.4i | %.4i | %.4i\n", i, fineTime_2[i], tot_2[i]);
	}
	printf("==================\n");	
	printf("\n==========================\n");
	printf("\n");

	return 1;
}



// Set the event to read verbose Ch = 1, 2, or 3 (both)
int EventTDC::SetVerbose(int ch){
	verbose = ch;
	return 1;
}




// Return the number of signals in one trigger Ch = 1, or 2
int EventTDC::size(int ch){
	int ret = 0;
	if(ch == 1){
		if(fineTime_1.size() == 0){
			ret = 0;
		}else if(tot_1.size() == 0){ 
			ret = 0;
		} else { ret = (int) fineTime_1.size();}
	} else if(ch == 2){
		if(fineTime_2.size() == 0){ 
			ret = 0;
		}else if(tot_2.size() == 0){ 
			ret = 0;
		} else { ret = (int) fineTime_2.size();}
	}
	return ret; 
}




// Search within the array for the search term - Ch 1 
int EventTDC::GetTDCTime1(){

	int count = 0;
	bin_data[0] = '\0';

	// Assign the binary data values using a bit mask
	int j = 0;
	for(int i=1; i >=0; i--){
		bin_data[(8 * i) + 0] = !!(channel_data[j] & 0b10000000);
		bin_data[(8 * i) + 1] = !!(channel_data[j] & 0b01000000);
		bin_data[(8 * i) + 2] = !!(channel_data[j] & 0b00100000);
		bin_data[(8 * i) + 3] = !!(channel_data[j] & 0b00010000);

		bin_data[(8 * i) + 4] = !!(channel_data[j] & 0b00001000);
		bin_data[(8 * i) + 5] = !!(channel_data[j] & 0b00000100);
		bin_data[(8 * i) + 6] = !!(channel_data[j] & 0b00000010);
		bin_data[(8 * i) + 7] = !!(channel_data[j] & 0b00000001);
		j = j+1;
	}

	// Search the binary data for the "search term" i.e. a change in value
	for(int i=15; i >= 4; i--){	
		if(bin_data[i] != searchTerm_1){
			count += 1;
		} else {
			delay_1 += count;

			if(searchTerm_1 == 1){
				if(Ch1_ready == 1){
					fineTime_1.push_back(delay_1);
					Ch1_ready = 0;
				}
				delay_1 = 0;
				count = 1;
				searchTerm_1 = 0;

			// When searching for the end of the signal
			} else if (searchTerm_1 == 0){
				if(Ch1_ready == 0){
					tot_1_candidate += delay_1;
				}
				delay_1 = 0;
				count = 1;
				searchTerm_1 = 1;
			}
		}
	}

	// If event is verbose then print the data
	if(verbose == 1 or verbose == 3){
		printf("Ch1 = [ ");
		for(int j = 0; j < 4; j++) printf("%u", bin_data[j]);
		printf(" ");
		for(int j = 4; j < 8; j++) printf("%u", bin_data[j]);
		printf(" ");
		for(int j = 8; j < 12; j++) printf("%u", bin_data[j]);
		printf(" ");
		for(int j = 12; j < 16; j++) printf("%u", bin_data[j]);
		printf("] :: %i :: %i\n", delay_1, (delay_1/12)); 
	}

	// Update the count value
	count = delay_1 + count;
	delay_1 = count;

	channel_data[0] = '\0';		
	return 1;
}



// Search within the array for the search term - Ch 2 
int EventTDC::GetTDCTime2(){

	int count = 0;
	bin_data[0] = '\0';

	// Assign the binary data values using a bit mask
	int j = 0;
	for(int i=1; i >=0; i--){
		bin_data[(8 * i) + 0] = !!(channel_data[j] & 0b10000000);
		bin_data[(8 * i) + 1] = !!(channel_data[j] & 0b01000000);
		bin_data[(8 * i) + 2] = !!(channel_data[j] & 0b00100000);
		bin_data[(8 * i) + 3] = !!(channel_data[j] & 0b00010000);

		bin_data[(8 * i) + 4] = !!(channel_data[j] & 0b00001000);
		bin_data[(8 * i) + 5] = !!(channel_data[j] & 0b00000100);
		bin_data[(8 * i) + 6] = !!(channel_data[j] & 0b00000010);
		bin_data[(8 * i) + 7] = !!(channel_data[j] & 0b00000001);
		j = j+1;
	}


	// Search the binary data for the "search term" i.e. a change in value
	for(int i=15; i >= 4; i--){	
		if(bin_data[i] != searchTerm_2){
			count += 1;
		} else {
			// If the data includes the search term save the value
			delay_2 += count;

			// When searching for the initial time
			if(searchTerm_2 == 1){
				if(Ch2_ready == 1){
					fineTime_2.push_back(delay_2);
					Ch2_ready = 0;
				}
				delay_2 = 0;
				count = 1;
				searchTerm_2 = 0;

			// When searching for the end of the signal
			} else if (searchTerm_2 == 0){
				if(Ch2_ready == 0){
					tot_2_candidate += delay_2;
				}
				delay_2 = 0;
				count = 1;
				searchTerm_2 = 1;
			}
		}
	}

	// If event is verbose then print the data
	if(verbose == 2 or verbose == 3){
		printf("Ch2 = [ ");
		for(int j = 0; j < 4; j++) printf("%u", bin_data[j]);
		printf(" ");
		for(int j = 4; j < 8; j++) printf("%u", bin_data[j]);
		printf(" ");
		for(int j = 8; j < 12; j++) printf("%u", bin_data[j]);
		printf(" ");
		for(int j = 12; j < 16; j++) printf("%u", bin_data[j]);
		printf("] :: %i :: %i\n", delay_2, (delay_2/12)); 
	}

	// If the data includes the search term save the value
	count = delay_2 + count;
	delay_2 = count;

	channel_data[0] = '\0';		
	return 1;
}






// Read an event from a file given
int EventTDC::Read(FILE *infile){ 
	
	int BytesRead = 0;
	
	// Read the payload length
	payloadLength = 0;
	BytesRead = fread(&payloadLength, sizeof(short int), 1, infile);	
	BytesRead -= 1;
	
	// Set starting search terms 
	searchTerm_1 = 1;
	delay_1 = 0; 
	searchTerm_2 = 1;
	delay_2 = 0; 
	Ch1_ready = 1;
	Ch2_ready = 1;
	noTriggers = 0;

	// Read the payload
	int LoopLen = payloadLength;
	for(int i = 0; i < LoopLen; i++){  


		//printf("i TDCsignalLen = %i\n", (i%TDCsignalLen));
		if(i%TDCsignalLen == 0){
			searchTerm_1 = 1;
			delay_1 = 0; 
			searchTerm_2 = 1;
			delay_2 = 0; 
			Ch1_ready = 1;
			Ch2_ready = 1;
			noTriggers += 1;
			tot_1_candidate = 0;
			tot_2_candidate = 0;
		}
		
		// Read 32 bits of data
		BytesRead += fread(&channel_data_A, sizeof(unsigned char), 2, infile);	
		BytesRead -= 2;
		BytesRead += fread(&channel_data, sizeof(unsigned char), 2, infile);	
		BytesRead -= 2;


		// Swap word order.
		GetTDCTime2();


		for(int bdc = 0; bdc < 12; bdc++){
			DataFrame2[(i%TDCsignalLen)*24 + bdc] += bin_data[(15 - bdc)];
		}

		memcpy(channel_data, channel_data_A, 2);
		GetTDCTime2();			
		for(int bdc = 0; bdc < 12; bdc++){
			DataFrame2[(i%TDCsignalLen)*24 + bdc + 12] += bin_data[(15 - bdc)];
		}

		channel_data_A[0] = '\0';

		// Read 32 bits of data
		BytesRead += fread(&channel_data_A, sizeof(unsigned char), 2, infile);	
		BytesRead -= 2;
		BytesRead += fread(&channel_data, sizeof(unsigned char), 2, infile);	
		BytesRead -= 2;

		GetTDCTime1();
		for(int bdc = 0; bdc < 12; bdc++){
			DataFrame1[(i%TDCsignalLen)*24 + bdc] += bin_data[(15 - bdc)];
		}
		memcpy(channel_data, channel_data_A, 2);
		GetTDCTime1();
		for(int bdc = 0; bdc < 12; bdc++){
			DataFrame1[(i%TDCsignalLen)*24 + bdc + 12] += bin_data[(15 - bdc)];
		}


		if(i%TDCsignalLen == TDCsignalLen-1){
			if(i > 0){
				id = 0;
  				BytesRead += fread(&id, sizeof(unsigned int), 1, infile);	
  				BytesRead -= 1;
				trigNo_2.push_back(id);
				BytesRead += fread(&id, sizeof(unsigned int), 1, infile);	
				BytesRead -= 1;
				trigNo_1.push_back(id);
				LoopLen -= 1;
			}
		}
		channel_data_A[0] = '\0';
	
		if(i%TDCsignalLen == TDCsignalLen-1){
			
			if(Ch1_ready == 0) tot_1.push_back(tot_1_candidate);
			if(Ch2_ready == 0) tot_2.push_back(tot_2_candidate);

			if(Ch1_ready == 1 and delay_1 > 0 and searchTerm_1 == 1){
				fineTime_1.push_back(-1);
				tot_1.push_back(-1);
			} else if(Ch1_ready == 0 and fineTime_1.size() != tot_1.size()){
				tot_1.push_back(delay_1);
			}

			if(Ch2_ready == 1 and delay_2 > 0 and searchTerm_2 == 1){
				fineTime_2.push_back(-1);
				tot_2.push_back(-1);
			} else if(Ch2_ready == 0 and fineTime_2.size() != tot_2.size()){
				tot_2.push_back(delay_2);
			}
		
		}

		// Set latch incase of runaway
		if(delay_1 > READ_LIMIT or delay_2 > READ_LIMIT){
			printf("READ_LIMIT reached!\n");
			printf("If this is too low change in file or investigate in verbose mode.\n");
			break;
		}

	}

	return BytesRead;
}






// Read a single file
std::vector<EventTDC> ReadChannel(char* filename, char* folderName, long int &file_loc){

	// Define basepath for reading in a folder
	char basepath[128];
	basepath[0] = '\0';
	strcat(basepath, folderName);
	strcat(basepath, "/DaqData/");
	strcat(basepath, filename);

	//============================================= Open File
	FILE *infile  = fopen(basepath, "rb"); 
	printf("Opened file %s\n", basepath);
	if(infile == NULL){
		printf("Can't open file '%s'!\n", filename);
	}	
	printf("Opened file %s\n", filename);

	//============================================= Read Header
	int BytesTotal = 0;
	int BytesRead = 0;
	
	HeaderTDC header = HeaderTDC(0); // HeaderTDC( Channel number ) 	
	BytesRead = header.Read(infile);
	header.Print();

	//============================================= Read Data
	std::vector<EventTDC> events;

	int ch_no;
	char hold[2], hold2[1] = {'\0'};
	int offset = strlen(filename);
	char hold3[1] = {'\0'};
	hold3[0] = filename[offset - 6]; 

	if(strcmp(hold3, "0123456789") == 1){
			hold[0] = filename[offset - 6];
			hold[1] = filename[offset - 5];
			ch_no = atoi((char*)hold);
	} else {
			hold2[0] = filename[offset - 5];
			ch_no = atoi(hold2);	
	}

	ch_no = ch_no * 2;
	// Iterate through events
	for(int i=0; i <= EVENT_LIMIT; i++){
		
		// Read events
		EventTDC event = EventTDC(i, ch_no); // EventTDC( id , ch_no) (int) 
		event.SetChID(ch_no);

		BytesRead = event.Read(infile);	

		if(BytesRead != 0){
			printf("Bytes finished at %i\n", i);
			break;
		}

		if(event.payloadLength == 0){
			event.noTriggers = 15;
			number_of_zero_len_payloads += 1;
		}

		events.push_back(event);
		
		// Set latch incase of runaway
		if(i > EVENT_LIMIT){
			printf("EVENT_LIMIT reached!\n");
			printf("If this is too low change in file or investigate in verbose mode.\n");
			break;
		}
	}
	
	printf("Read file\n");
	
	file_loc = ftell(infile);
	fclose(infile);
	
	return events;
}

// Read a single file
long int HistoGUI::ReReadChannel(char* filename, char* folderName, long int seek_to){

	// Define basepath for reading in a folder
	char basepath[128];
	basepath[0] = '\0';
	strcat(basepath, folderName);
	strcat(basepath, "/DaqData/");
	strcat(basepath, filename);

	//============================================= Open File
	FILE *infile  = fopen(basepath, "rb"); 
//	printf("Opened file %s\n", basepath);
	if(infile == NULL){
		printf("Can't open file '%s'!\n", filename);
	}	
//	printf("Opened file %s\n", filename);
	fseek(infile, seek_to, SEEK_SET);

	//============================================= Read Header - NO NEED
	int BytesTotal = 0;
	int BytesRead = 0;
	

	//============================================= Read Data

	int ch_no;
	char hold[2], hold2[1] = {'\0'};
	int offset = strlen(filename);
	char hold3[1] = {'\0'};
	hold3[0] = filename[offset - 6]; 

	if(strcmp(hold3, "0123456789") == 1){
			hold[0] = filename[offset - 6];
			hold[1] = filename[offset - 5];
			ch_no = atoi((char*)hold);
	} else {
			hold2[0] = filename[offset - 5];
			ch_no = atoi(hold2);	
	}

	ch_no = ch_no * 2;
	// Iterate through events
	for(int i=0; i <= EVENT_LIMIT; i++){
		
		// Read events
		EventTDC event = EventTDC(i, ch_no); // EventTDC( id , ch_no) (int) 
		event.SetChID(ch_no);

		BytesRead = event.Read(infile);	

		if(BytesRead != 0){
			printf("Bytes finished at %i\n", i);
			break;
		}

		if(event.payloadLength == 0){
			event.noTriggers = 15;
			number_of_zero_len_payloads += 1;
		}


		for(int j=0; j < event.size(1); j++){
			if(event.fineTime_1[j] > 0){
					y[event.fineTime_1[j]] += 1;
					content_2[event.fineTime_1[j]][event.tot_1[j]] += 1;
					content_3[div(64 - pixelmap[ch_no], 8).rem -4 + 1][div(64 - pixelmap[ch_no], 8).quot + 1] += 1;
			}
			if(event.fineTime_2[j] > 0){
					y[event.fineTime_2[j]] += 1;
					content_2[event.fineTime_2[j]][event.tot_2[j]] += 1;
					content_3[div(64 - pixelmap[ch_no +1], 8).rem -4 + 1][div(64 - pixelmap[ch_no +1], 8).quot + 1] += 1;
			}
		}
		
		// Set latch incase of runaway
		if(i > EVENT_LIMIT){
//			printf("EVENT_LIMIT reached!\n");
//			printf("If this is too low change in file or investigate in verbose mode.\n");
			break;
		}
	}
	
	printf("Read file\n");
	
	long int file_loc = ftell(infile);
	fclose(infile);
	
	return file_loc;
}


// ======= ==== ========
// ======== MAIN ========
// ======== ==== ========


#ifndef __CINT__ 

int main(int argc,char **argv){


	// READ FOLDER =============================================================
	DIR *folder;
	char fpath[128];
	sprintf(fpath, "%s/DaqData", argv[1]);
	folder = opendir(fpath);

	// Error message
	if(folder == NULL){
	    printf("Unable to read '%s'!\n",argv[1]);
	    return(1);
	} else {
	    printf("Reading from '%s'.\n", argv[1]);
	}

	struct dirent *entry;
	int folderStart = telldir(folder);

	int files = 0;
	while( (entry=readdir(folder))) {
		if(strcmp(entry->d_name, ".") == 0) continue;
		if(strcmp(entry->d_name, "..") == 0) continue;
		files++;
	}

	// Return to the top of the file
	seekdir(folder,folderStart);
		

	// MAP =============================================================
    std::map<int, int> pixelmap = {
        //ChNo  PIXEL
		// H - Type
		{ 0 , 52 ,},
		{ 1 , 57 ,},
		{ 2 , 60 ,},
		{ 3 , 58 ,},
		{ 4 , 59 ,},
		{ 5 , 50 ,}, 
		{ 6 , 49 ,}, 
		{ 7 , 51 ,},
		{ 8 , 36 ,},  
		{ 9 , 41 ,},
		{ 10, 44 ,}, 
		{ 11, 42 ,},
		{ 12, 43 ,},
		{ 13, 33 ,},
		{ 14, 27 ,},
		{ 15, 34 ,},
		{ 16, 35 ,},
		{ 17, 18 ,},
		{ 18, 19 ,},
		{ 19, 25 ,},
		{ 20, 26 ,},
		{ 21, 28 ,},
		{ 22, 9  ,},
		{ 23, 10 ,},
		{ 24, 20 ,},
		{ 25, 1  ,},
		{ 26, 3  ,},
		{ 27, 12 ,},
		{ 28, 11 ,},
		{ 29, 4  ,},
		{ 30, 2  ,},
		{ 31, 17 ,}
	};

	double RisingEdgesArray[24 * TDCsignalLen][32] = {0};
	std::vector<double> RisingEdgesVect(24 * TDCsignalLen);
	double HitMapArray[8][8] = {0};
	double ChannelOccupancy[32] = {0};


	// LOOP TRHOUGH FILES ========================================================

	// Counters ---------
	int currentFile = 0;
	int noRisingEdges = 0;
	std::vector<std::vector<int>> All_channels; //(32);

	printf("READING ======================================\n\n");
	long int file_loc[16] = {0};

	

	HistoGUI gui;
	for(int cn=0; cn < 24 * TDCsignalLen; cn++){
		gui.x.push_back(cn);
		gui.y.push_back(0);

		gui.x_2.push_back(cn);
		gui.y_2.push_back(cn);
		std::vector<double> temp;
		for(int ca=0; ca < 24 * TDCsignalLen; ca++){
			temp.push_back(0);	
		}
		gui.content_2.push_back(temp);
	}
	for(int cn=0; cn < 8; cn++){
		gui.x_3.push_back(cn);
		gui.y_3.push_back(cn);

		std::vector<double> temp;
		for(int ca=0; ca < 8; ca++){
			temp.push_back(0);	
		}
		gui.content_3.push_back(temp);
	}
	char all_filenames[16][64];

	// Loop -------------
	while( (entry=readdir(folder))) {
		if(strcmp(entry->d_name, ".") == 0) continue;
		if(strcmp(entry->d_name, "..") == 0) continue;

		// GET CHANNEL NAME =======================================
		int ch_no = 0;
		char hold[3], hold2[2];
		hold[2] = '\0';
		hold2[1] = '\0';
		int offset = strlen(entry->d_name);
		char hold3[1] = {'\0'};
		hold3[0] = entry->d_name[offset - 6]; 

		if(strcmp(hold3, "0123456789") == 1){
				hold[0] = entry->d_name[offset - 6];
				hold[1] = entry->d_name[offset - 5];
				ch_no = atoi((char*)hold);
			//	printf("Channel is double digit: %c %c\n", entry->d_name[offset - 6], entry->d_name[offset - 5]);
		} else {
			//	printf("Channel is single digit: %c\n", entry->d_name[offset - 5]);
				hold2[0] = entry->d_name[offset - 5];
				ch_no = atoi(hold2);	
			//	printf("hold2 = %s -> %i\n", hold2, ch_no);
		}

		ch_no = ch_no * 2;

		// Read File ------------
		std::vector <EventTDC> events = ReadChannel(entry->d_name, argv[1], file_loc[(ch_no / 2)]);
		sprintf(all_filenames[(ch_no / 2)], "%s", entry->d_name);
		std::vector <int> channel;
		printf("%s : Loc = %i\n", entry->d_name, file_loc[(ch_no / 2)]);

		// Alert if empty -------
		if(events.size() == 0){
			printf("FILE %s is empty!\n", entry->d_name);
			continue;
		}

		// ITERATE TRHOUGH DATA =======================================
		unsigned long long int count_ch1 = 0;
		unsigned long long int count_ch2 = 0;	
		int payload_sum = 0;
		int noTriggers_channel = 0;
		int noTriggers_channel_2 = 0;

		printf("Ch No = %i\n", ch_no);
		printf("No Payloads = %lu\n", events.size());

		// Loop --------------
		for(int i = 0; i < events.size(); i++){
			payload_sum += events[i].payloadLength;
			channel.push_back(events[i].payloadLength);
			noTriggers_channel += events[i].noTriggers;
			noTriggers_channel_2 += events[i].size(1);
			if(events.size() == 0) break;
			//printf("%i : %i\n", i, events[i].payloadLength);

			// Loop through triggers --------------
			for(int j=0; j < events[i].size(1); j++){

				// If ToA <= 0, set to -1 & fill tree
				if(events[i].fineTime_1[j] <= 0){
					continue;
				} else {

					count_ch1 += 1;
        			HitMapArray[div(64 - pixelmap[ch_no], 8).rem -4 + 1][div(64 - pixelmap[ch_no], 8).quot + 1] += 1;//, (1.0/noTriggers));
					//if(id)
					RisingEdgesArray[events[i].fineTime_1[j]][ch_no] += 1;
					gui.y[events[i].fineTime_1[j]] += 1;
					gui.content_2[events[i].fineTime_1[j]][events[i].tot_1[j]] += 1;
					gui.content_3[div(64 - pixelmap[ch_no], 8).rem -4 + 1][div(64 - pixelmap[ch_no], 8).quot + 1] += 1;
					// Increment counters ---------------
					noRisingEdges += 1;
	
				}
			}
		}


		ch_no = ch_no+1;
		printf("Ch No = %i\n", ch_no);
		for(int i = 0; i < events.size(); i++){
			for(int j=0; j < events[i].size(2); j++){

				if(events[i].fineTime_2[j] <= 0){
					continue;
				}else{

					count_ch2++;
					noRisingEdges += 1;
					
        			HitMapArray[div(64 - pixelmap[ch_no], 8).rem -4][div(64 - pixelmap[ch_no], 8).quot] += 1;//, (1.0/noTriggers));
					RisingEdgesArray[events[i].fineTime_2[j]][ch_no] += 1;
					gui.y[events[i].fineTime_2[j]] += 1;
					gui.content_2[events[i].fineTime_2[j]][events[i].tot_2[j]] += 1;
					gui.content_3[div(64 - pixelmap[ch_no], 8).rem -4 + 1][div(64 - pixelmap[ch_no], 8).quot + 1] += 1;
				}

			}
		}	

		printf("Payload sum   = %i\n", payload_sum);
		printf("Number of zero length payloads = %i\n", number_of_zero_len_payloads);
		printf("\nDATA:\n");
		printf("Ch %.2i = %llu / %i = %.3f%% \n", ch_no-1, count_ch1, noTriggers_channel, 100. * (double)count_ch1/noTriggers_channel );
		ChannelOccupancy[ch_no-1] = (double)count_ch1 / noTriggers_channel;
		printf("Ch %.2i = %llu / %i = %.3f%% \n", ch_no, count_ch2, noTriggers_channel, 100. * (double)count_ch2/noTriggers_channel );
		ChannelOccupancy[ch_no] = (double)count_ch2 / noTriggers_channel;
		printf("No Triggers   = %i\n", noTriggers_channel);
//		printf("No Triggers 2 = %i\n", noTriggers_channel_2);
		printf("Final Trig count = %lu\n", events[(events.size()-1)].id);
//		printf("First Trig count = %lu\n", events[0].trigNo_1[0]);
	//	printf("Trig diff        = %lu\n", events[(events.size()-1)].id - events[0].trigNo_1[0]);
		printf("\n");
		number_of_zero_len_payloads = 0;
		All_channels.push_back(channel);	
		printf("\n---------------------------------------------- \n\n");
	}	

	gui.Init();	
	for(int kl=0;kl<16; kl++) printf("%s\n", all_filenames[kl]);
	gui.TDCLoop(file_loc, all_filenames, argv[1]);
	gui.Close();

	printf("Total number of rising edges = %i\n", noRisingEdges);

	int hm_sum = 0;
	for(int k = 0; k < 8; k++){
		for(int m = 0; m < 8; m++){
			hm_sum += HitMapArray[k][m];
		}
	}

	printf("\n\n\n");
	printf("GRAPHS =======================================\n\n");


	printf("Rough Hitmap - Not mapped!\n");
	printf("--------------------------\n\n");
	for(int k = 0; k < 8; k++){
		printf(" | ");
		for(int m = 0; m < 8; m++){
			//printf("%.2f  |  ", HitMapArray[k][m] / noRisingEdges);
			if(HitMapArray[m][k] / hm_sum > 0.05){
			printf("X | ");
			} else if(HitMapArray[m][k] / hm_sum > 0.01){
			printf("+ | ");
			} else {
			printf("  | ");
			}
		}
		printf("\n");
	}

	
	printf("\n\n\n");
	
	printf("Rough Timing Alignment - Not optimised!\n");
	printf("---------------------------------------\n\n");

	int re_sum[32] = {0};
	double means[32] = {0};
	for(int k = 0; k < 32; k++){
		for(int m = 0; m < 24 * TDCsignalLen; m++){
			re_sum[k] += RisingEdgesArray[m][k];
		}
		means[k] = re_sum[k] / (24 * TDCsignalLen);
	}

	for(int c = 0; c < 32; c++){
		printf("Ch %.2i = %.3f%%  : ", c, 100. * ChannelOccupancy[c]);
		for(int k = 0; k < 24; k++){
			int junk_count = 0;
			for(int m = 0; m < TDCsignalLen; m++){
				junk_count += RisingEdgesArray[m + (k * TDCsignalLen)][c];
			}
				
			if((double) junk_count / re_sum[c] > 0.20){
				printf("X");
			} else if((double) junk_count / re_sum[c] > 0.10){
				printf("+");
			} else {
				printf("-");
			}	
		}
		printf("\n");
	}

	printf("\n\n\n");




//	for (int i = 0; i < All_channels.size(); i++){
//		int pl_sum_2 = 0;
//		printf("Channel = %i\n",i);
//		for (int j = 0; j < All_channels[i].size(); j++){
//			printf("%.3i ", All_channels[i][j]);
//			pl_sum_2 += All_channels[i][j];
//			if((j+1)%8 == 0) printf("\n");
//		}
//		printf("payload sum = %i\n", pl_sum_2);
//		printf("\n\n");
//	}




	return 1;

}

#endif


/*
int main(int argc,char **argv){

	HistoGUI gui;

	std::vector<double> x;
	std::vector<double> y;

	for(int i = 0; i < 25; i++){
		x.push_back( i );
		y.push_back( pow(i, 0.5) );
	}

	gui.SetData(x,y);
	gui.Init();
	gui.Loop();
	gui.Close();

	return 1;

}

*/









