#include "global.h"

//To compile : g++ AD.cxx -o {Input Executable Name} -lX11
// VERSION : 1.1 (Overlay for Angular Fit Display) 
// Created with thanks to Daniel Folds Holt - University of Cambridge Cavendish Lab




using namespace std;


class HistoGUIad{


	public:

		HistoGUIad(){}

		Display * disp;
		Window    wind;
		XEvent    evt;
		int       screen;

		XColor xcolour;
		XColor xcolour_red;
		XColor xcolour_veryred;
		Colormap cmap;


		std::vector<double> x;
		std::vector<double> y;
		std::vector<double> y_errors;

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


		double A0;
		double A2E;
		double A4E; 
		int SetFit(double a, double b, double c){
			A0 = a; 
			A2E = b;
			A4E = c;
			fit_exists = true;
			return 1;
		}


		double old_xl, old_xh, old_yl, old_yh;
		double old_mouse_x, old_mouse_y;

		int Init();
		int SetData(std::vector<double> a, std::vector<double>b);
		int SetErrors(std::vector<double> a);
		int Loop();	
		void Close(){ printf("Closing now!\n"); XCloseDisplay(disp); }
		int DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
		int Draw_Fit(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
		int DrawCrosshairs(int mouse_x, int mouse_y);
		int Zoom(int mouse_x, int mouse_y);

		double legval(double theta);
		bool fit_exists = false;


//		int Legendre_fit(std::vector<double> d1, std::vector<double> d2, double 

};

int HistoGUIad::SetData(std::vector<double> a, std::vector<double>b){
	for(int i=0; i < a.size(); i++){
		x.push_back(a[i]);
		y.push_back(b[i]);
	}
	return a.size();
}

int HistoGUIad::SetErrors(std::vector<double> a){
	for(int i=0; i < a.size(); i++){
		y_errors.push_back(a[i]);
	}
	return a.size();
}
int HistoGUIad::Init(){

	disp = XOpenDisplay(NULL);
	if (disp == NULL) {
		fprintf(stderr, "Cannot open display\n");
		exit(1);
	}

	//fit_exists = false;

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



	return 1;

}

int HistoGUIad::DrawCrosshairs(int mouse_x, int mouse_y){

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

int HistoGUIad::Zoom(int mouse_x, int mouse_y){

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
double HistoGUIad::legval(double theta){
	
	double lg;
	
	double	aaa = A0/A0; 
	double 	aab  = (A2E/A0)*(1.5 * pow(cos(theta),2) - .5);
    double  aac  = (A4E/A0)*(35./8. * pow(cos(theta),4) - 30./8. * pow(cos(theta),2)  +  3./8. );
	
	
	lg = aaa + aab + aac;


	return lg;
}

int HistoGUIad::Draw_Fit(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win){


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

		width_scale = (max_x - min_x) / (0.8 * width);
		x_offset = 0.5 * ((max_x - min_x) / 0.8) - 0.5 * (min_x + max_x);

		height_scale = -1. * (max_y - min_y) / (0.8 * height);
		y_offset = -0.5 * ((max_y - min_y) / 0.8) + -0.5 * (min_y + max_y);

		XSetForeground(disp, DefaultGC(disp,screen), xcolour_red.pixel);

		double x_wid ;
		double y_wid ;

		for(int i=0; i < (int) width; i += 2){
			double x_val   = i * width_scale - x_offset;
			double y_val = legval(x_val);
			
			y_wid  = (y_val + y_offset) / height_scale;
			XFillRectangle(disp, wind, DefaultGC(disp, screen), i -2, y_wid -2, 4, 4);
		}
	//	XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid2 -2, y_wid2 -2, 4, 4);

	} else {
		//double x_low_win, double y_low_win, double x_hi_win, double y_hi_win

		width_scale = (x_hi_win - x_low_win) / width;
		x_offset = -1.0 * x_low_win;

		height_scale = (y_hi_win - y_low_win) / height;
		y_offset = -1.0 * y_low_win;


		XSetForeground(disp, DefaultGC(disp,screen), xcolour_red.pixel);

		double x_wid ;
		double y_wid ;

		for(int i=0; i < (int) width; i += 2){
			double x_val   = i * width_scale - x_offset;
			double y_val = legval(x_val);
			
			y_wid  = (y_val + y_offset) / height_scale;
			XFillRectangle(disp, wind, DefaultGC(disp, screen), i -2, y_wid -2, 4, 4);
		}
	}

	old_xl = x_low_win;
	old_yl = y_low_win;
	old_xh = x_hi_win;
	old_yh = y_hi_win;



	return 1; 
}

int HistoGUIad::DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win){
		
	if(fit_exists){
		Draw_Fit(x_low_win,y_low_win,x_hi_win,y_hi_win); 

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
	//	double x_wid2;
	//	double y_wid2;
		double y_errors_wid;

		for(int i=0; i < x.size() - 1; i++){
			x_wid  = (x[i] + x_offset) / width_scale;
			y_wid  = (y[i] + y_offset) / height_scale;
			y_errors_wid = (y_errors[i] / height_scale);
	//		x_wid2 = (x[i + 1] + x_offset) / width_scale;
	//		y_wid2 = (y[i + 1] + y_offset) / height_scale;
			//printf("(%f, %f), (%f,%f)\n", x_wid,y_wid,x_wid2,y_wid2);
	//		XDrawLine(disp, wind, DefaultGC(disp, screen), x_wid, y_wid, x_wid2, y_wid2);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -3, y_wid -3, 6, 6);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -1, y_wid + y_errors_wid, 2, y_errors_wid);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -1, y_wid - y_errors_wid, 2, y_errors_wid);
//error bar testing 		
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -4, y_wid + y_errors_wid, 8, 2);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -4, y_wid - y_errors_wid, 8, 2);
		}
	//	XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid2 -2, y_wid2 -2, 4, 4);

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
	//	double x_wid2;
	//	double y_wid2;

		double y_errors_wid;
		
		for(int i=0; i < x.size() - 1; i++){
			x_wid  = (x[i] + x_offset) / width_scale;
			y_wid  = (y[i] + y_offset) / height_scale;
	
			y_errors_wid = (y_errors[i] / height_scale);

		//	x_wid2 = (x[i + 1] + x_offset) / width_scale;
		//	y_wid2 = (y[i + 1] + y_offset) / height_scale;
			//	printf("(%f, %f), (%f,%f)\n", x_wid,y_wid,x_wid2,y_wid2);
		//	XDrawLine(disp, wind, DefaultGC(disp, screen), x_wid, y_wid, x_wid2, y_wid2);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -3, y_wid -3, 6, 6);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -1, y_wid , 2, y_errors_wid);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -1, y_wid - y_errors_wid, 2, y_errors_wid);
//testig for error bars.		
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -4, y_wid + y_errors_wid, 8, 2);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid -4, y_wid - y_errors_wid, 8, 2);
		}
	//	XFillRectangle(disp, wind, DefaultGC(disp, screen), x_wid2 -2, y_wid2 -2, 4, 4);


	}

	old_xl = x_low_win;
	old_yl = y_low_win;
	old_xh = x_hi_win;
	old_yh = y_hi_win;

	return 1;
}


int HistoGUIad::Loop(){

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
				DrawCrosshairs(mouse_x, mouse_y);
				DrawData(old_xl, old_yl, old_xh, old_yh);
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

			DrawCrosshairs(old_mouse_x, old_mouse_y);
			DrawCrosshairs(mouse_x, mouse_y);
			DrawData(old_xl, old_yl, old_xh, old_yh);

		} else if (evt.type == MotionNotify and MousePressed2){

			int mouse_x = evt.xmotion.x;
			int mouse_y = evt.xmotion.y;

			XSetForeground(disp, DefaultGC(disp,screen), xcolour_veryred.pixel);
			XFillRectangle(disp, wind, DefaultGC(disp, screen), mouse_x -2, mouse_y -2, 4, 4);

		} else if (evt.type == KeyPress){

			printf("Key pressed = %x\n", evt.xkey.keycode);
			if(evt.xkey.keycode == 0x41 or evt.xkey.keycode == 0x39){
				XClearWindow(disp, wind);
				DrawData(-1,-1,-1,-1);

			} else {
				break;	
			}
		}
	}

	return 1;
}
