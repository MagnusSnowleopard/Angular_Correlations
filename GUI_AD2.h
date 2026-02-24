#ifndef GUI_AD_H
#define GUI_AD_H

//#############################################################################
//
// VERSION : 1.3 (Improved Angular Fit Display GUI)
//
// Purpose:
//   Display angular distributions from CSV input with error bars and fit overlay.
//   Improved with axis labels, plot margins, larger fonts, safer buffers,
//   consistent zoom, and cleaner styling.
//
// Styling requested:
//   - All lines/text black
//   - Fit line red
//   - Bigger fonts (~2x)
//
//#############################################################################

#include "global.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <X11/keysym.h>

class HistoGUIad {
public:
    HistoGUIad();

    // Core API (compatible with your current calls)
    int Init();
    int SetData(const std::vector<double>& a, const std::vector<double>& b);
    int SetErrors(const std::vector<double>& a);
    int SetFit(double a, double b, double c);
    int Loop();
    void Close();

    int DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
    int Draw_Fit(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
    int DrawCrosshairs(int mouse_x, int mouse_y);
    int Zoom(int mouse_x, int mouse_y);

    double legval(double theta);

    // Optional polish/configuration
    void SetLabels(const std::string& xlab, const std::string& ylab, const std::string& title = "");
    void SetMargins(int left, int right, int top, int bottom);

    bool fit_exists;

private:
    // X11
    Display* disp;
    Window   wind;
    XEvent   evt;
    int      screen;

    XColor xcolour_red;      // fit (red)
    XColor xcolour_veryred;  // paint mode (kept)
    Colormap cmap;

    // Fonts (bigger)
    XFontStruct* font_big;
    XFontStruct* font_small;

    // Data
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> y_errors;

    // Current data window bounds (data coords)
    double max_x;
    double min_x;
    double max_y;
    double min_y;

    // Mouse positions / converted coords
    double pos_x;
    double pos_y;
    double new_pos_x;
    double new_pos_y;

    // Window geometry
    unsigned int width;
    unsigned int height;

    // Data<->pixel transform:
    // data_x = px * width_scale  - x_offset
    // data_y = py * height_scale - y_offset
    double width_scale;
    double x_offset;
    double height_scale;
    double y_offset;

    // Fit params
    double A0;
    double A2E;
    double A4E;
    double Ierr;

    // Last view (for redraw/zoom)
    double old_xl, old_xh, old_yl, old_yh;
    double old_mouse_x, old_mouse_y;

    // Labels / styling
    std::string plot_title;
    std::string x_axis_label;
    std::string y_axis_label;

    int margin_left;
    int margin_right;
    int margin_top;
    int margin_bottom;

    // Helpers
    int  UpdateGeometry();
    void ComputeDefaultBounds(double& xl, double& yl, double& xh, double& yh);
    void SetTransform(double xl, double yl, double xh, double yh);
    void DrawAxesAndLabels();
    bool HasUsableData() const;

    // Plot rectangle helpers
    int PlotLeft()   const;
    int PlotRight()  const;
    int PlotTop()    const;
    int PlotBottom() const;
    int PlotW()      const;
    int PlotH()      const;

    // Transform helpers
    double PixToDataX(int px) const;
    double PixToDataY(int py) const;
    int    DataToPixX(double xv) const;
    int    DataToPixY(double yv) const;

    bool InPlot(int px, int py) const;
};

//=============================
// Constructor
//=============================
inline HistoGUIad::HistoGUIad()
    : fit_exists(false),
      disp(NULL),
      wind(0),
      screen(0),
      cmap(0),
      font_big(NULL),
      font_small(NULL),
      max_x(1.0), min_x(0.0), max_y(1.0), min_y(0.0),
      pos_x(0.0), pos_y(0.0), new_pos_x(0.0), new_pos_y(0.0),
      width(900), height(600),
      width_scale(1.0), x_offset(0.0), height_scale(-1.0), y_offset(0.0),
      A0(1.0), A2E(0.0), A4E(0.0), Ierr(0.0),
      old_xl(-1.0), old_xh(-1.0), old_yl(-1.0), old_yh(-1.0),
      old_mouse_x(0.0), old_mouse_y(0.0),
      plot_title("Angular Distribution + Fit"),
      x_axis_label("Theta (rad)"),
      y_axis_label("W(theta) / A0"),
      margin_left(95), margin_right(20), margin_top(35), margin_bottom(70)
{}

//=============================
// Public config methods
//=============================
inline void HistoGUIad::SetLabels(const std::string& xlab,
                                  const std::string& ylab,
                                  const std::string& title)
{
    x_axis_label = xlab;
    y_axis_label = ylab;
    if (!title.empty()) plot_title = title;
}

inline void HistoGUIad::SetMargins(int left, int right, int top, int bottom)
{
    margin_left   = std::max(40, left);
    margin_right  = std::max(5,  right);
    margin_top    = std::max(20, top);
    margin_bottom = std::max(35, bottom);
}

//=============================
// Data setters
//=============================
inline int HistoGUIad::SetData(const std::vector<double>& a, const std::vector<double>& b)
{
    x.clear();
    y.clear();

    size_t n = std::min(a.size(), b.size());
    x.reserve(n);
    y.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        x.push_back(a[i]);
        y.push_back(b[i]);
    }

    if (y_errors.size() != x.size()) {
        y_errors.assign(x.size(), 0.0);
    }

    return static_cast<int>(n);
}

inline int HistoGUIad::SetErrors(const std::vector<double>& a)
{
    y_errors.assign(x.size(), 0.0);
    size_t n = std::min(a.size(), x.size());
    for (size_t i = 0; i < n; ++i) {
        y_errors[i] = a[i];
    }
    return static_cast<int>(n);
}

inline int HistoGUIad::SetFit(double a, double b, double c)
{
    A0 = a;
    A2E = b;
    A4E = c;
    fit_exists = true;
    return 1;
}

//=============================
// X11 init/close
//=============================
inline int HistoGUIad::Init()
{
    disp = XOpenDisplay(NULL);
    if (disp == NULL) {
        std::fprintf(stderr, "Cannot open display\n");
        std::exit(1);
    }

    screen = DefaultScreen(disp);

    wind = XCreateSimpleWindow(
        disp,
        RootWindow(disp, screen),
        10, 10,
        900, 600,
        1,
        BlackPixel(disp, screen),
        WhitePixel(disp, screen)
    );

    XGrabPointer(disp, wind, False,
                 ButtonPressMask | ButtonReleaseMask | PointerMotionMask,
                 GrabModeAsync, GrabModeAsync,
                 None, None, CurrentTime);

    XSelectInput(disp, wind,
                 ExposureMask | KeyPressMask |
                 ButtonPressMask | ButtonReleaseMask |
                 PointerMotionMask | StructureNotifyMask);

    XMapWindow(disp, wind);

    // Red fit color
    cmap = DefaultColormap(disp, screen);

    xcolour_red.red   = 50000;
    xcolour_red.green = 0;
    xcolour_red.blue  = 0;
    xcolour_red.flags = DoRed | DoGreen | DoBlue;
    XAllocColor(disp, cmap, &xcolour_red);

    // paint mode color (kept)
    xcolour_veryred.red   = 65000;
    xcolour_veryred.green = 12000;
    xcolour_veryred.blue  = 12000;
    xcolour_veryred.flags = DoRed | DoGreen | DoBlue;
    XAllocColor(disp, cmap, &xcolour_veryred);

    // Load larger fonts (~2x)
    // Fallback chain varies by system
    font_big = XLoadQueryFont(disp, "-misc-fixed-bold-r-normal--20-200-75-75-c-100-iso8859-1");
    if (!font_big) font_big = XLoadQueryFont(disp, "10x20");
    if (!font_big) font_big = XLoadQueryFont(disp, "9x15bold");
    if (!font_big) font_big = XQueryFont(disp, XGContextFromGC(DefaultGC(disp, screen)));

    font_small = XLoadQueryFont(disp, "9x15");
    if (!font_small) font_small = font_big;

    if (font_big) {
        XSetFont(disp, DefaultGC(disp, screen), font_big->fid);
    }

    // Slightly thicker lines help readability
    XSetLineAttributes(disp, DefaultGC(disp, screen), 2, LineSolid, CapButt, JoinMiter);

    return 1;
}

inline void HistoGUIad::Close()
{
    std::printf("Closing now!\n");

    if (disp) {
        XFontStruct* fb = font_big;
        XFontStruct* fs = font_small;

        if (fb) XFreeFont(disp, fb);
        if (fs && fs != fb) XFreeFont(disp, fs);

        font_big = NULL;
        font_small = NULL;

        XCloseDisplay(disp);
        disp = NULL;
    }
}

//=============================
// Geometry / transform helpers
//=============================
inline int HistoGUIad::UpdateGeometry()
{
    int j1, j2;
    unsigned int j3, j4;
    Window root_return;

    XGetGeometry(disp, wind, &root_return, &j1, &j2, &width, &height, &j3, &j4);
    return 1;
}

inline int HistoGUIad::PlotLeft() const   { return margin_left; }
inline int HistoGUIad::PlotRight() const  { return (int)width - margin_right; }
inline int HistoGUIad::PlotTop() const    { return margin_top; }
inline int HistoGUIad::PlotBottom() const { return (int)height - margin_bottom; }

inline int HistoGUIad::PlotW() const
{
    int w = PlotRight() - PlotLeft();
    return (w > 1) ? w : 1;
}

inline int HistoGUIad::PlotH() const
{
    int h = PlotBottom() - PlotTop();
    return (h > 1) ? h : 1;
}

inline double HistoGUIad::PixToDataX(int px) const
{
    return px * width_scale - x_offset;
}

inline double HistoGUIad::PixToDataY(int py) const
{
    return py * height_scale - y_offset;
}

inline int HistoGUIad::DataToPixX(double xv) const
{
    return (int)((xv + x_offset) / width_scale);
}

inline int HistoGUIad::DataToPixY(double yv) const
{
    return (int)((yv + y_offset) / height_scale);
}

inline bool HistoGUIad::InPlot(int px, int py) const
{
    return (px >= PlotLeft() && px <= PlotRight() &&
            py >= PlotTop()  && py <= PlotBottom());
}

inline bool HistoGUIad::HasUsableData() const
{
    return (!x.empty() && x.size() == y.size());
}

inline void HistoGUIad::SetTransform(double xl, double yl, double xh, double yh)
{
    if (std::fabs(xh - xl) < 1e-12) xh = xl + 1.0;
    if (std::fabs(yh - yl) < 1e-12) yh = yl + 1.0;

    // x mapping: px=PlotLeft -> xl, px=PlotRight -> xh
    width_scale = (xh - xl) / (double)PlotW();
    x_offset    = PlotLeft() * width_scale - xl;

    // y mapping: py=PlotTop -> yh, py=PlotBottom -> yl (math-style y up)
    height_scale = -(yh - yl) / (double)PlotH();
    y_offset     = PlotTop() * height_scale - yh;

    old_xl = xl; old_xh = xh;
    old_yl = yl; old_yh = yh;
}

inline void HistoGUIad::ComputeDefaultBounds(double& xl, double& yl, double& xh, double& yh)
{
    if (!HasUsableData()) {
        xl = 0.0;
        xh = 3.141592653589793;
        yl = 0.0;
        yh = 1.0;
        return;
    }

    double x_min = x[0];
    double x_max = x[0];
    double y_min = y[0];
    double y_max = y[0];

    for (size_t i = 0; i < x.size(); ++i) {
        x_min = std::min(x_min, x[i]);
        x_max = std::max(x_max, x[i]);

        double yi = y[i];
        double ei = (i < y_errors.size()) ? y_errors[i] : 0.0;

        y_min = std::min(y_min, yi - ei);
        y_max = std::max(y_max, yi + ei);
    }

    // Include fit envelope in autoscale if fit exists
    if (fit_exists) {
        int nsamp = 200;
        double xs0 = x_min;
        double xs1 = x_max;
        if (std::fabs(xs1 - xs0) < 1e-12) xs1 = xs0 + 1.0;

        for (int i = 0; i <= nsamp; ++i) {
            double t = (double)i / (double)nsamp;
            double xv = xs0 + t * (xs1 - xs0);
            double fv = legval(xv);
            y_min = std::min(y_min, fv);
            y_max = std::max(y_max, fv);
        }
    }

    double x_pad = 0.05 * (x_max - x_min);
    double y_pad = 0.10 * (y_max - y_min);

    if (std::fabs(x_pad) < 1e-12) x_pad = 0.5;
    if (std::fabs(y_pad) < 1e-12) y_pad = 0.5;

    xl = x_min - x_pad;
    xh = x_max + x_pad;
    yl = y_min - y_pad;
    yh = y_max + y_pad;
}

//=============================
// Fit function
//=============================
inline double HistoGUIad::legval(double theta)
{
    // Uses your normalized form (scaled by A0 in coefficients)
    double norm = (std::fabs(A0) > 1e-14) ? A0 : 1.0;

    double c  = std::cos(theta);
    double c2 = c * c;
    double c4 = c2 * c2;

    double aaa = A0 / norm;
    double aab = (A2E / norm) * (1.5 * c2 - 0.5);
    double aac = (A4E / norm) * (35.0/8.0 * c4 - 30.0/8.0 * c2 + 3.0/8.0);

    return (aaa + aab + aac);
}

//=============================
// Axes / labels
//=============================
inline void HistoGUIad::DrawAxesAndLabels()
{
    // ALL non-fit lines/text black (requested)
    XSetForeground(disp, DefaultGC(disp, screen), BlackPixel(disp, screen));

    // Plot rectangle
    XDrawRectangle(disp, wind, DefaultGC(disp, screen),
                   PlotLeft(), PlotTop(), PlotW(), PlotH());

    // Optional x=0 and y=0 axes if visible
    int axis_x = DataToPixX(0.0);
    int axis_y = DataToPixY(0.0);

    if (axis_x >= PlotLeft() && axis_x <= PlotRight()) {
        XDrawLine(disp, wind, DefaultGC(disp, screen),
                  axis_x, PlotTop(), axis_x, PlotBottom());
    }
    if (axis_y >= PlotTop() && axis_y <= PlotBottom()) {
        XDrawLine(disp, wind, DefaultGC(disp, screen),
                  PlotLeft(), axis_y, PlotRight(), axis_y);
    }

    // Title + axis labels in bigger font
    if (font_big) XSetFont(disp, DefaultGC(disp, screen), font_big->fid);

    // Title
    XDrawString(disp, wind, DefaultGC(disp, screen),
                PlotLeft() + 8, 22,
                plot_title.c_str(), (int)plot_title.size());

    // X-axis label (bottom center-ish)
    int xlab_x = PlotLeft() + PlotW()/2 - (int)(x_axis_label.size() * 5);
    if (xlab_x < PlotLeft()) xlab_x = PlotLeft();
    XDrawString(disp, wind, DefaultGC(disp, screen),
                xlab_x, (int)height - 18,
                x_axis_label.c_str(), (int)x_axis_label.size());

    // Y-axis label (horizontal due to Xlib text limitations)
    XDrawString(disp, wind, DefaultGC(disp, screen),
                8, 22,
                y_axis_label.c_str(), (int)y_axis_label.size());

    // Tick labels (smaller but still larger than default)
    if (font_small) XSetFont(disp, DefaultGC(disp, screen), font_small->fid);

    const int nx = 8;
    const int ny = 8;
    char axis_val[64];

    for (int i = 0; i <= nx; ++i) {
        int px = PlotLeft() + (i * PlotW()) / nx;
        double xv = PixToDataX(px);

        XDrawLine(disp, wind, DefaultGC(disp, screen), px, PlotBottom(), px, PlotBottom() + 6);

        std::snprintf(axis_val, sizeof(axis_val), "%.2f", xv);
        XDrawString(disp, wind, DefaultGC(disp, screen),
                    px - 16, PlotBottom() + 22, axis_val, (int)std::strlen(axis_val));
    }

    for (int i = 0; i <= ny; ++i) {
        int py = PlotTop() + (i * PlotH()) / ny;
        double yv = PixToDataY(py);

        XDrawLine(disp, wind, DefaultGC(disp, screen), PlotLeft() - 6, py, PlotLeft(), py);

        std::snprintf(axis_val, sizeof(axis_val), "%.2f", yv);
        XDrawString(disp, wind, DefaultGC(disp, screen),
                    8, py + 5, axis_val, (int)std::strlen(axis_val));
    }

    // Restore big font as default for any later labels/crosshair text
    if (font_big) XSetFont(disp, DefaultGC(disp, screen), font_big->fid);
}

//=============================
// Crosshairs
//=============================
inline int HistoGUIad::DrawCrosshairs(int mouse_x, int mouse_y)
{
    UpdateGeometry();

    // Clamp to plot region so crosshairs are visually clean
    int cx = std::max(PlotLeft(), std::min(mouse_x, PlotRight()));
    int cy = std::max(PlotTop(),  std::min(mouse_y, PlotBottom()));

    // Crosshairs black (requested: all non-fit lines black)
    XSetForeground(disp, DefaultGC(disp, screen), BlackPixel(disp, screen));
    if (font_big) XSetFont(disp, DefaultGC(disp, screen), font_big->fid);

    XDrawLine(disp, wind, DefaultGC(disp, screen), cx, PlotTop(), cx, PlotBottom());
    XDrawLine(disp, wind, DefaultGC(disp, screen), PlotLeft(), cy, PlotRight(), cy);

    pos_x = PixToDataX(cx);
    pos_y = PixToDataY(cy);

    char coord[128];
    std::snprintf(coord, sizeof(coord), "(%.4f, %.4f)", pos_x, pos_y);

    int tx = cx + 10;
    int ty = cy - 10;
    if (tx > (int)width - 220) tx = cx - 180;
    if (ty < 20) ty = cy + 24;

    XDrawString(disp, wind, DefaultGC(disp, screen), tx, ty, coord, (int)std::strlen(coord));

    return 1;
}

//=============================
// Zoom
//=============================
inline int HistoGUIad::Zoom(int mouse_x, int mouse_y)
{
    // Clamp drag release to plot region
    int mx = std::max(PlotLeft(), std::min(mouse_x, PlotRight()));
    int my = std::max(PlotTop(),  std::min(mouse_y, PlotBottom()));

    new_pos_x = PixToDataX(mx);
    new_pos_y = PixToDataY(my);

    int ox = std::max(PlotLeft(), std::min((int)old_mouse_x, PlotRight()));
    int oy = std::max(PlotTop(),  std::min((int)old_mouse_y, PlotBottom()));

    pos_x = PixToDataX(ox);
    pos_y = PixToDataY(oy);

    // Ignore tiny drags
    if (std::abs(mx - ox) < 12 || std::abs(my - oy) < 12) {
        return 1;
    }

    double x_low_win = std::min(pos_x, new_pos_x);
    double x_hi_win  = std::max(pos_x, new_pos_x);
    double y_low_win = std::min(pos_y, new_pos_y);
    double y_hi_win  = std::max(pos_y, new_pos_y);

    XClearWindow(disp, wind);
    DrawData(x_low_win, y_low_win, x_hi_win, y_hi_win);
    return 1;
}

//=============================
// Fit drawing (red only)
//=============================
inline int HistoGUIad::Draw_Fit(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win)
{
    if (!fit_exists) return 1;

    UpdateGeometry();

    // If called independently, establish a transform
    if (x_low_win == -1.0 && y_low_win == -1.0 && x_hi_win == -1.0 && y_hi_win == -1.0) {
        double xl, yl, xh, yh;
        ComputeDefaultBounds(xl, yl, xh, yh);
        SetTransform(xl, yl, xh, yh);
    }

    // Fit line red (requested)
    XSetForeground(disp, DefaultGC(disp, screen), xcolour_red.pixel);

    bool have_prev = false;
    int prev_px = 0, prev_py = 0;

    for (int px = PlotLeft(); px <= PlotRight(); ++px) {
        double x_val = PixToDataX(px);
        double y_val = legval(x_val);
        int py = DataToPixY(y_val);

        if (!InPlot(px, py)) {
            have_prev = false;
            continue;
        }

        if (have_prev) {
            XDrawLine(disp, wind, DefaultGC(disp, screen), prev_px, prev_py, px, py);
        } else {
            XFillRectangle(disp, wind, DefaultGC(disp, screen), px - 1, py - 1, 3, 3);
        }

        prev_px = px;
        prev_py = py;
        have_prev = true;
    }

    return 1;
}

//=============================
// Main data drawing
//=============================
inline int HistoGUIad::DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win)
{
    UpdateGeometry();

    if (!HasUsableData()) {
        XSetForeground(disp, DefaultGC(disp, screen), BlackPixel(disp, screen));
        if (font_big) XSetFont(disp, DefaultGC(disp, screen), font_big->fid);
        const char* msg = "No data loaded for HistoGUIad";
        XDrawString(disp, wind, DefaultGC(disp, screen), 20, 40, msg, (int)std::strlen(msg));
        return 1;
    }

    // Resolve bounds (autoscale or zoom window)
    double xl, yl, xh, yh;
    if (x_low_win == -1.0 && y_low_win == -1.0 && x_hi_win == -1.0 && y_hi_win == -1.0) {
        ComputeDefaultBounds(xl, yl, xh, yh);
    } else {
        xl = x_low_win;
        yl = y_low_win;
        xh = x_hi_win;
        yh = y_hi_win;
    }

    min_x = xl; min_y = yl;
    max_x = xh; max_y = yh;

    SetTransform(xl, yl, xh, yh);

    // Draw axes/box/labels/ticks in black
    DrawAxesAndLabels();

    // Fit overlay first (red)
    if (fit_exists) {
        Draw_Fit(xl, yl, xh, yh);
    }

    // Data + error bars in black
    XSetForeground(disp, DefaultGC(disp, screen), BlackPixel(disp, screen));

    for (size_t i = 0; i < x.size(); ++i) {
        int px = DataToPixX(x[i]);
        int py = DataToPixY(y[i]);

        if (px < PlotLeft() || px > PlotRight()) continue;

        // Keep compatibility with your original error scaling behavior
        double e = 0.0;
        if (i < y_errors.size()) {
            if (std::fabs(A0) > 1e-12) e = y_errors[i] / A0;
            else e = y_errors[i];
        }

        int e_pix = (int)std::fabs(e / height_scale);

        if (py < PlotTop() - 20 || py > PlotBottom() + 20) continue;

        // Vertical error bar
        XDrawLine(disp, wind, DefaultGC(disp, screen),
                  px, py - e_pix, px, py + e_pix);

        // Caps
        XDrawLine(disp, wind, DefaultGC(disp, screen),
                  px - 4, py - e_pix, px + 4, py - e_pix);
        XDrawLine(disp, wind, DefaultGC(disp, screen),
                  px - 4, py + e_pix, px + 4, py + e_pix);

        // Point marker
        XFillRectangle(disp, wind, DefaultGC(disp, screen), px - 3, py - 3, 6, 6);
    }

    return 1;
}

//=============================
// Event loop
//=============================
inline int HistoGUIad::Loop()
{
    bool MousePressed  = false; // left drag = crosshair + zoom
    bool MousePressed2 = false; // right drag = paint mode (kept)

    while (1) {
        XNextEvent(disp, &evt);

        if (evt.type == Expose) {
            DrawData(-1, -1, -1, -1);

        } else if (evt.type == ConfigureNotify) {
            // Window resized
            XClearWindow(disp, wind);
            if (old_xl == -1.0 && old_yl == -1.0 && old_xh == -1.0 && old_yh == -1.0)
                DrawData(-1, -1, -1, -1);
            else
                DrawData(old_xl, old_yl, old_xh, old_yh);

        } else if (evt.type == ButtonPress) {
            int mouse_x = evt.xbutton.x;
            int mouse_y = evt.xbutton.y;
            old_mouse_x = mouse_x;
            old_mouse_y = mouse_y;

            if (evt.xbutton.button == Button1) {
                // Left click/drag for crosshair + zoom
                XClearWindow(disp, wind);

                if (old_xl == -1.0 && old_yl == -1.0 && old_xh == -1.0 && old_yh == -1.0)
                    DrawData(-1, -1, -1, -1);
                else
                    DrawData(old_xl, old_yl, old_xh, old_yh);

                // Draw crosshair AFTER data so it remains visible
                DrawCrosshairs(mouse_x, mouse_y);

                MousePressed = true;
                std::printf("(%.6f, %.6f)\n", pos_x, pos_y);
            }

            if (evt.xbutton.button == Button3) {
                MousePressed2 = true;
            }

            if (evt.xbutton.button == Button2) {
                // Middle button reset view
                XClearWindow(disp, wind);
                DrawData(-1, -1, -1, -1);
            }

        } else if (evt.type == ButtonRelease) {
            int mouse_x = evt.xbutton.x;
            int mouse_y = evt.xbutton.y;

            if (evt.xbutton.button == Button1) {
                Zoom(mouse_x, mouse_y);
                MousePressed = false;
            }

            if (evt.xbutton.button == Button3) {
                MousePressed2 = false;
            }

        } else if (evt.type == MotionNotify && MousePressed) {
            int mouse_x = evt.xmotion.x;
            int mouse_y = evt.xmotion.y;

            XClearWindow(disp, wind);

            if (old_xl == -1.0 && old_yl == -1.0 && old_xh == -1.0 && old_yh == -1.0)
                DrawData(-1, -1, -1, -1);
            else
                DrawData(old_xl, old_yl, old_xh, old_yh);

            DrawCrosshairs((int)old_mouse_x, (int)old_mouse_y);
            DrawCrosshairs(mouse_x, mouse_y);

        } else if (evt.type == MotionNotify && MousePressed2) {
            // Paint mode kept, in black or red? (using very-red marker)
            int mouse_x = evt.xmotion.x;
            int mouse_y = evt.xmotion.y;

            XSetForeground(disp, DefaultGC(disp, screen), xcolour_veryred.pixel);
            XFillRectangle(disp, wind, DefaultGC(disp, screen), mouse_x - 2, mouse_y - 2, 4, 4);

        } else if (evt.type == KeyPress) {
            KeySym ks = XLookupKeysym(&evt.xkey, 0);

            // Reset keys: r, R, space
            if (ks == XK_r || ks == XK_R || ks == XK_space) {
                XClearWindow(disp, wind);
                DrawData(-1, -1, -1, -1);
                continue;
            }

            // Quit keys: q, Q, Esc
            if (ks == XK_q || ks == XK_Q || ks == XK_Escape) {
                break;
            }

            // Any other key -> ignore
        }
    }

    return 1;
}

#endif // GUI_AD_H
