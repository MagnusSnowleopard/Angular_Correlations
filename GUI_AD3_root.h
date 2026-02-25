#ifndef GUI_AD_ROOT_H
#define GUI_AD_ROOT_H

//#############################################################################
//
// ROOT-based Angular Distribution GUI (replacement for X11 version)
//
// Features:
//   - ROOT GUI (TGMainFrame) + embedded TCanvas
//   - Top controls: j1, j2, E_gamma (keV), file path, Browse/Load buttons
//   - Option buttons 1 / 2 / 3 on top toolbar
//   - TGraphErrors for angular data + error bars
//   - Optional fit overlay (Legendre A0 + A2 P2 + A4 P4), drawn in red
//   - ROOT native mouse zoom/pan on canvas (no manual X11 crosshair loop)
//
// Notes:
//   - Option buttons 1/2/3 are wired and tracked (fOptionMode), but you can
//     define their physics-specific behavior later.
//   - File parser accepts comma, space, tab (mixed delimiters supported).
//   - Expects at least 3 columns per data row: angle, y, ey
//
// Compile example (standalone):
//   g++ -std=c++17 your_main.cpp `root-config --cflags --glibs` -lGui -lGuiHtml -o gui_ad
//
// Typical use:
//   auto* gui = new HistoGUIadRoot(gClient->GetRoot(), 1200, 800);
//   gui->Init();
//   gui->Loop();   // optional if standalone app and you want blocking behavior
//
//#############################################################################

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// ROOT GUI / graphics
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGLayout.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>

//------------------------------------------------------------------------------
// HistoGUIadRoot
//------------------------------------------------------------------------------
class HistoGUIadRoot : public TGMainFrame {
public:
    HistoGUIadRoot(const TGWindow* p = gClient->GetRoot(),
                   UInt_t w = 1200,
                   UInt_t h = 800);

    virtual ~HistoGUIadRoot();

    // Core API (kept similar to your original)
    int Init();
    int SetData(const std::vector<double>& a, const std::vector<double>& b);
    int SetErrors(const std::vector<double>& a);
    int SetFit(double a, double b, double c);
    int Loop();   // optional blocking loop for standalone app
    void Close();

    // Kept for compatibility; ROOT now handles interactive zoom/pan in canvas
    int DrawData(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
    int Draw_Fit(double x_low_win, double y_low_win, double x_hi_win, double y_hi_win);
    int DrawCrosshairs(int /*mouse_x*/, int /*mouse_y*/) { return 1; } // ROOT-native interaction
    int Zoom(int /*mouse_x*/, int /*mouse_y*/) { return 1; }            // ROOT-native interaction

    double legval(double theta);

    // Optional polish/configuration
    void SetLabels(const std::string& xlab, const std::string& ylab, const std::string& title = "");
    void SetMargins(double left, double right, double top, double bottom);

    // GUI parameter accessors
    double GetJ1() const;
    double GetJ2() const;
    double GetGammaKeV() const;
    std::string GetFileName() const;

    void SetJ1(double v);
    void SetJ2(double v);
    void SetGammaKeV(double v);
    void SetFileName(const std::string& s);

    // Load file from GUI file entry (or explicit path)
    bool LoadDataFromFile(const std::string& filename);
    bool LoadFromGuiFileAndDraw();

    // ROOT GUI command handler
    Bool_t ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t parm2) override;

    // Window close handler
    void CloseWindow() override;

    bool fit_exists;

private:
    enum EWidgetIds {
        kID_Browse = 1001,
        kID_Load   = 1002,
        kID_Reset  = 1003,
        kID_Quit   = 1004,
        kID_Opt1   = 1101,
        kID_Opt2   = 1102,
        kID_Opt3   = 1103,
        kID_Redraw = 1201
    };

    // GUI builders/helpers
    void BuildGUI();
    void UpdateModeButtonStates();
    void UpdateTitleText();
    void ResetView();
    void RefreshPlot();

    // Parsing helpers
    static std::vector<std::string> TokenizeFlexible(const std::string& line);
    static bool TryParse3(const std::vector<std::string>& row, double& a, double& b, double& c);

    // Plot helpers
    bool HasUsableData() const;
    void ComputeDefaultBounds(double& xl, double& yl, double& xh, double& yh) const;
    void BuildOrUpdateGraphs();
    void ApplyPadMargins();
    void StyleFrameHistogram(TH1* hframe);

private:
    // GUI widgets
    TGVerticalFrame*        fMainV;
    TGHorizontalFrame*      fTopRow1;
    TGHorizontalFrame*      fTopRow2;
    TRootEmbeddedCanvas*    fECanvas;

    TGLabel*                fLblJ1;
    TGLabel*                fLblJ2;
    TGLabel*                fLblEg;
    TGLabel*                fLblFile;

    TGNumberEntry*          fEntJ1;
    TGNumberEntry*          fEntJ2;
    TGNumberEntry*          fEntEg;

    TGTextEntry*            fEntFile;

    TGTextButton*           fBtnBrowse;
    TGTextButton*           fBtnLoad;
    TGTextButton*           fBtnReset;
    TGTextButton*           fBtnQuit;

    TGTextButton*           fBtnOpt1;
    TGTextButton*           fBtnOpt2;
    TGTextButton*           fBtnOpt3;
    TGTextButton*           fBtnRedraw;

    // Data
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> y_errors;

    // Optional raw content if you want to preserve original behavior
    std::vector< std::vector<std::string> > content;

    // Plot objects (owned by class)
    TGraphErrors*           fGraph;
    TGraph*                 fFitGraph;
    TLegend*                fLegend;

    // Frame pointer from DrawFrame (owned by pad/canvas, not deleted here)
    TH1*                    fFrameHist;

    // Fit params
    double A0;
    double A2E;
    double A4E;
    double Ierr;

    // Last drawn view (for explicit redraw compatibility)
    double old_xl, old_xh, old_yl, old_yh;

    // Data bounds cache (not strictly required but kept for compatibility)
    double max_x;
    double min_x;
    double max_y;
    double min_y;

    // Labels / styling
    std::string plot_title;
    std::string x_axis_label;
    std::string y_axis_label;

    // ROOT pad margins (fractional)
    double margin_left;
    double margin_right;
    double margin_top;
    double margin_bottom;

    // Misc GUI state
    int fOptionMode;  // 1,2,3
    bool fIsMapped;
};

//=============================
// Constructor
//=============================
inline HistoGUIadRoot::HistoGUIadRoot(const TGWindow* p, UInt_t w, UInt_t h)
    : TGMainFrame(p, w, h),
      fit_exists(false),
      fMainV(nullptr),
      fTopRow1(nullptr),
      fTopRow2(nullptr),
      fECanvas(nullptr),
      fLblJ1(nullptr),
      fLblJ2(nullptr),
      fLblEg(nullptr),
      fLblFile(nullptr),
      fEntJ1(nullptr),
      fEntJ2(nullptr),
      fEntEg(nullptr),
      fEntFile(nullptr),
      fBtnBrowse(nullptr),
      fBtnLoad(nullptr),
      fBtnReset(nullptr),
      fBtnQuit(nullptr),
      fBtnOpt1(nullptr),
      fBtnOpt2(nullptr),
      fBtnOpt3(nullptr),
      fBtnRedraw(nullptr),
      fGraph(nullptr),
      fFitGraph(nullptr),
      fLegend(nullptr),
      fFrameHist(nullptr),
      A0(1.0), A2E(0.0), A4E(0.0), Ierr(0.0),
      old_xl(-1.0), old_xh(-1.0), old_yl(-1.0), old_yh(-1.0),
      max_x(1.0), min_x(0.0), max_y(1.0), min_y(0.0),
      plot_title("Angular Distribution + Fit"),
      x_axis_label("Theta (rad)"),
      y_axis_label("W(#theta) / A_{0}"),
      margin_left(0.12), margin_right(0.04), margin_top(0.08), margin_bottom(0.12),
      fOptionMode(1),
      fIsMapped(false)
{
    SetCleanup(kDeepCleanup); // ROOT will recursively clean child widgets
    BuildGUI();
    UpdateModeButtonStates();
    UpdateTitleText();
}

//=============================
// Destructor
//=============================
inline HistoGUIadRoot::~HistoGUIadRoot()
{
    delete fGraph;    fGraph = nullptr;
    delete fFitGraph; fFitGraph = nullptr;
    delete fLegend;   fLegend = nullptr;

    // Child TG widgets are cleaned by kDeepCleanup on TGMainFrame cleanup
    Cleanup();
}

//=============================
// GUI builders/helpers
//=============================
inline void HistoGUIadRoot::BuildGUI()
{
    // Main vertical frame
    fMainV = new TGVerticalFrame(this);

    // ---------------------------
    // Row 1: physics inputs + file
    // ---------------------------
    fTopRow1 = new TGHorizontalFrame(fMainV);

    fLblJ1 = new TGLabel(fTopRow1, "j1:");
    fEntJ1 = new TGNumberEntry(fTopRow1, 0.0, 8, -1,
                               TGNumberFormat::kNESReal,
                               TGNumberFormat::kNEANonNegative,
                               TGNumberFormat::kNELLimitMinMax,
                               0.0, 100.0);

    fLblJ2 = new TGLabel(fTopRow1, "j2:");
    fEntJ2 = new TGNumberEntry(fTopRow1, 0.0, 8, -1,
                               TGNumberFormat::kNESReal,
                               TGNumberFormat::kNEANonNegative,
                               TGNumberFormat::kNELLimitMinMax,
                               0.0, 100.0);

    fLblEg = new TGLabel(fTopRow1, "E#gamma (keV):");
    fEntEg = new TGNumberEntry(fTopRow1, 0.0, 10, -1,
                               TGNumberFormat::kNESReal,
                               TGNumberFormat::kNEANonNegative,
                               TGNumberFormat::kNELLimitMinMax,
                               0.0, 1.0e9);

    fLblFile = new TGLabel(fTopRow1, "File:");
    fEntFile = new TGTextEntry(fTopRow1, "");

    fBtnBrowse = new TGTextButton(fTopRow1, "&Browse", kID_Browse);
    fBtnLoad   = new TGTextButton(fTopRow1, "&Load / Plot", kID_Load);

    fBtnBrowse->Associate(this);
    fBtnLoad->Associate(this);

    fTopRow1->AddFrame(fLblJ1,    new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 6, 4, 4, 4));
    fTopRow1->AddFrame(fEntJ1,    new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 8, 4, 4));

    fTopRow1->AddFrame(fLblJ2,    new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 6, 4, 4, 4));
    fTopRow1->AddFrame(fEntJ2,    new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 8, 4, 4));

    fTopRow1->AddFrame(fLblEg,    new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 6, 4, 4, 4));
    fTopRow1->AddFrame(fEntEg,    new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 12, 4, 4));

    fTopRow1->AddFrame(fLblFile,  new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 6, 4, 4, 4));
    fTopRow1->AddFrame(fEntFile,  new TGLayoutHints(kLHintsCenterY | kLHintsExpandX, 0, 8, 4, 4));

    fTopRow1->AddFrame(fBtnBrowse,new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 4, 4, 4));
    fTopRow1->AddFrame(fBtnLoad,  new TGLayoutHints(kLHintsCenterY | kLHintsRight, 2, 6, 4, 4));

    // ---------------------------
    // Row 2: option buttons + actions
    // ---------------------------
    fTopRow2 = new TGHorizontalFrame(fMainV);

    fBtnOpt1  = new TGTextButton(fTopRow2, "1", kID_Opt1);
    fBtnOpt2  = new TGTextButton(fTopRow2, "2", kID_Opt2);
    fBtnOpt3  = new TGTextButton(fTopRow2, "3", kID_Opt3);
    fBtnRedraw= new TGTextButton(fTopRow2, "Redraw", kID_Redraw);
    fBtnReset = new TGTextButton(fTopRow2, "Reset Zoom", kID_Reset);
    fBtnQuit  = new TGTextButton(fTopRow2, "Quit", kID_Quit);

    fBtnOpt1->Associate(this);
    fBtnOpt2->Associate(this);
    fBtnOpt3->Associate(this);
    fBtnRedraw->Associate(this);
    fBtnReset->Associate(this);
    fBtnQuit->Associate(this);

    fTopRow2->AddFrame(fBtnOpt1,   new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 2, 4, 4));
    fTopRow2->AddFrame(fBtnOpt2,   new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 4, 4));
    fTopRow2->AddFrame(fBtnOpt3,   new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 12, 4, 4));

    fTopRow2->AddFrame(fBtnRedraw, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 4, 4, 4));
    fTopRow2->AddFrame(fBtnReset,  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 4, 4, 4));

    fTopRow2->AddFrame(new TGLabel(fTopRow2, " "), new TGLayoutHints(kLHintsExpandX)); // spacer

    fTopRow2->AddFrame(fBtnQuit,   new TGLayoutHints(kLHintsRight | kLHintsCenterY, 4, 6, 4, 4));

    // Embedded ROOT canvas
    fECanvas = new TRootEmbeddedCanvas("ad_canvas", fMainV, 1000, 650);

    // Add rows/canvas to main frame
    fMainV->AddFrame(fTopRow1, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    fMainV->AddFrame(fTopRow2, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
    fMainV->AddFrame(fECanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 6, 6, 4, 6));

    AddFrame(fMainV, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // Window setup
    SetWindowName("Angular Distribution GUI (ROOT)");
    MapSubwindows();
    Resize(GetDefaultSize());

    // Canvas defaults
    if (fECanvas && fECanvas->GetCanvas()) {
        TCanvas* c = fECanvas->GetCanvas();
        c->SetFillColor(0);
        c->SetBorderMode(0);
        c->SetFrameBorderMode(0);
        c->SetGrid(0,0);
    }
}

inline void HistoGUIadRoot::UpdateModeButtonStates()
{
    // Simple visual cue via button labels (works across ROOT versions reliably)
    fBtnOpt1->SetText(fOptionMode == 1 ? "[1]" : "1");
    fBtnOpt2->SetText(fOptionMode == 2 ? "[2]" : "2");
    fBtnOpt3->SetText(fOptionMode == 3 ? "[3]" : "3");
}

inline void HistoGUIadRoot::UpdateTitleText()
{
    std::ostringstream os;
    os << "Angular Distribution GUI (ROOT)   |   "
       << "Mode=" << fOptionMode
       << "   j1=" << GetJ1()
       << "   j2=" << GetJ2()
       << "   E#gamma=" << GetGammaKeV() << " keV";
    SetWindowName(os.str().c_str());
}

inline void HistoGUIadRoot::ApplyPadMargins()
{
    if (!fECanvas || !fECanvas->GetCanvas()) return;
    TCanvas* c = fECanvas->GetCanvas();
    c->SetLeftMargin(margin_left);
    c->SetRightMargin(margin_right);
    c->SetTopMargin(margin_top);
    c->SetBottomMargin(margin_bottom);
}

inline void HistoGUIadRoot::StyleFrameHistogram(TH1* hframe)
{
    if (!hframe) return;

    hframe->SetTitle(plot_title.c_str());
    hframe->GetXaxis()->SetTitle(x_axis_label.c_str());
    hframe->GetYaxis()->SetTitle(y_axis_label.c_str());

    // Bigger fonts (~2x feel vs ROOT defaults)
    hframe->GetXaxis()->SetTitleSize(0.050);
    hframe->GetYaxis()->SetTitleSize(0.050);
    hframe->GetXaxis()->SetLabelSize(0.040);
    hframe->GetYaxis()->SetLabelSize(0.040);
    hframe->GetXaxis()->SetTitleOffset(1.05);
    hframe->GetYaxis()->SetTitleOffset(1.10);

    hframe->SetLineColor(kBlack);
    hframe->SetStats(0);
}

inline void HistoGUIadRoot::ResetView()
{
    old_xl = old_xh = old_yl = old_yh = -1.0;
    DrawData(-1, -1, -1, -1);
}

inline void HistoGUIadRoot::RefreshPlot()
{
    if (HasUsableData()) {
        DrawData(-1, -1, -1, -1);
    } else {
        // Draw empty frame with labels so GUI looks alive before loading
        if (!fECanvas || !fECanvas->GetCanvas()) return;
        TCanvas* c = fECanvas->GetCanvas();
        c->cd();
        c->Clear();
        ApplyPadMargins();
        fFrameHist = c->DrawFrame(0.0, 0.0, TMath::Pi(), 1.0);
        StyleFrameHistogram(fFrameHist);
        c->Modified();
        c->Update();
    }
}

//=============================
// Public config methods
//=============================
inline void HistoGUIadRoot::SetLabels(const std::string& xlab,
                                      const std::string& ylab,
                                      const std::string& title)
{
    x_axis_label = xlab;
    y_axis_label = ylab;
    if (!title.empty()) plot_title = title;
    RefreshPlot();
}

inline void HistoGUIadRoot::SetMargins(double left, double right, double top, double bottom)
{
    margin_left   = std::max(0.02, std::min(0.40, left));
    margin_right  = std::max(0.01, std::min(0.25, right));
    margin_top    = std::max(0.02, std::min(0.25, top));
    margin_bottom = std::max(0.05, std::min(0.35, bottom));
    RefreshPlot();
}

//=============================
// GUI parameter accessors
//=============================
inline double HistoGUIadRoot::GetJ1() const
{
    return fEntJ1 ? fEntJ1->GetNumber() : 0.0;
}

inline double HistoGUIadRoot::GetJ2() const
{
    return fEntJ2 ? fEntJ2->GetNumber() : 0.0;
}

inline double HistoGUIadRoot::GetGammaKeV() const
{
    return fEntEg ? fEntEg->GetNumber() : 0.0;
}

inline std::string HistoGUIadRoot::GetFileName() const
{
    return (fEntFile && fEntFile->GetText()) ? std::string(fEntFile->GetText()) : std::string();
}

inline void HistoGUIadRoot::SetJ1(double v)
{
    if (fEntJ1) fEntJ1->SetNumber(v);
    UpdateTitleText();
}

inline void HistoGUIadRoot::SetJ2(double v)
{
    if (fEntJ2) fEntJ2->SetNumber(v);
    UpdateTitleText();
}

inline void HistoGUIadRoot::SetGammaKeV(double v)
{
    if (fEntEg) fEntEg->SetNumber(v);
    UpdateTitleText();
}

inline void HistoGUIadRoot::SetFileName(const std::string& s)
{
    if (fEntFile) fEntFile->SetText(s.c_str());
}

//=============================
// Data setters
//=============================
inline int HistoGUIadRoot::SetData(const std::vector<double>& a, const std::vector<double>& b)
{
    x.clear();
    y.clear();

    const size_t n = std::min(a.size(), b.size());
    x.reserve(n);
    y.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        x.push_back(a[i]);
        y.push_back(b[i]);
    }

    if (y_errors.size() != x.size()) {
        y_errors.assign(x.size(), 0.0);
    }

    BuildOrUpdateGraphs();
    return static_cast<int>(n);
}

inline int HistoGUIadRoot::SetErrors(const std::vector<double>& a)
{
    y_errors.assign(x.size(), 0.0);
    const size_t n = std::min(a.size(), x.size());
    for (size_t i = 0; i < n; ++i) {
        y_errors[i] = a[i];
    }

    BuildOrUpdateGraphs();
    return static_cast<int>(n);
}

inline int HistoGUIadRoot::SetFit(double a, double b, double c)
{
    A0 = a;
    A2E = b;
    A4E = c;
    fit_exists = true;

    if (HasUsableData()) DrawData(-1, -1, -1, -1);
    return 1;
}

//=============================
// Init / Close / Loop
//=============================
inline int HistoGUIadRoot::Init()
{
    if (!fIsMapped) {
        MapWindow();
        fIsMapped = true;
    }
    RefreshPlot();
    return 1;
}

inline void HistoGUIadRoot::Close()
{
    // Standard ROOT GUI close path (assumes heap allocation, which is typical)
    SendCloseMessage();
}

inline void HistoGUIadRoot::CloseWindow()
{
    // Called when window manager close button is clicked
    DeleteWindow();
}

inline int HistoGUIadRoot::Loop()
{
    // Optional blocking loop for standalone usage.
    // If you're using ROOT interactively, you can skip calling this.
    if (gApplication) {
        gApplication->Run(kTRUE); // return when app is terminated
    }
    return 1;
}

//=============================
// Flexible parsing helpers
//=============================
inline std::vector<std::string> HistoGUIadRoot::TokenizeFlexible(const std::string& line)
{
    // Accept commas + whitespace (spaces/tabs), mixed in any combination
    std::string s = line;
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == ',') s[i] = ' ';
        // tabs/newlines are already whitespace for operator>>
    }

    std::vector<std::string> row;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok) row.push_back(tok);
    return row;
}

inline bool HistoGUIadRoot::TryParse3(const std::vector<std::string>& row, double& a, double& b, double& c)
{
    if (row.size() < 3) return false;

    char* e1 = nullptr;
    char* e2 = nullptr;
    char* e3 = nullptr;

    const double v1 = std::strtod(row[0].c_str(), &e1);
    const double v2 = std::strtod(row[1].c_str(), &e2);
    const double v3 = std::strtod(row[2].c_str(), &e3);

    if (e1 == row[0].c_str() || e2 == row[1].c_str() || e3 == row[2].c_str())
        return false;

    a = v1; b = v2; c = v3;
    return true;
}

//=============================
// File loading
//=============================
inline bool HistoGUIadRoot::LoadDataFromFile(const std::string& filename)
{
    if (filename.empty()) return false;

    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        std::fprintf(stderr, "Could not open file: %s\n", filename.c_str());
        return false;
    }

    content.clear();
    x.clear();
    y.clear();
    y_errors.clear();

    std::string line;
    while (std::getline(file, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Optional comment support
        std::string trimmed = line;
        trimmed.erase(trimmed.begin(),
                      std::find_if(trimmed.begin(), trimmed.end(),
                                   [](unsigned char ch){ return !std::isspace(ch); }));
        if (trimmed.empty()) continue;
        if (trimmed[0] == '#') continue;

        std::vector<std::string> row = TokenizeFlexible(line);
        if (row.empty()) continue;

        content.push_back(row);

        double aa = 0.0, yy = 0.0, ee = 0.0;
        if (!TryParse3(row, aa, yy, ee)) {
            // Skip header/non-numeric line safely
            continue;
        }

        x.push_back(aa);
        y.push_back(yy);
        y_errors.push_back(ee);
    }

    BuildOrUpdateGraphs();

    if (HasUsableData()) {
        DrawData(-1, -1, -1, -1);
        return true;
    }

    std::fprintf(stderr, "File loaded but no valid numeric rows found (need 3 columns).\n");
    RefreshPlot();
    return false;
}

inline bool HistoGUIadRoot::LoadFromGuiFileAndDraw()
{
    return LoadDataFromFile(GetFileName());
}

//=============================
// Plot data helpers
//=============================
inline bool HistoGUIadRoot::HasUsableData() const
{
    return (!x.empty() && x.size() == y.size());
}

inline void HistoGUIadRoot::BuildOrUpdateGraphs()
{
    delete fGraph;    fGraph = nullptr;
    delete fFitGraph; fFitGraph = nullptr;
    delete fLegend;   fLegend = nullptr;

    if (!HasUsableData()) return;

    const int n = static_cast<int>(x.size());
    fGraph = new TGraphErrors(n);

    for (int i = 0; i < n; ++i) {
        const double ey = (i < (int)y_errors.size()) ? y_errors[i] : 0.0;
        fGraph->SetPoint(i, x[(size_t)i], y[(size_t)i]);
        fGraph->SetPointError(i, 0.0, ey);
    }

    // Styling requested:
    //   - all non-fit lines/text black
    //   - fit line red
    fGraph->SetLineColor(kBlack);
    fGraph->SetMarkerColor(kBlack);
    fGraph->SetMarkerStyle(20);
    fGraph->SetMarkerSize(1.0);
    fGraph->SetLineWidth(2);
}

inline void HistoGUIadRoot::ComputeDefaultBounds(double& xl, double& yl, double& xh, double& yh) const
{
    if (!HasUsableData()) {
        xl = 0.0;
        xh = TMath::Pi();
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

        const double ei = (i < y_errors.size()) ? y_errors[i] : 0.0;
        y_min = std::min(y_min, y[i] - ei);
        y_max = std::max(y_max, y[i] + ei);
    }

    // Include fit envelope in autoscale if fit exists
    if (fit_exists) {
        const int nsamp = 500;
        double xs0 = x_min;
        double xs1 = x_max;
        if (std::fabs(xs1 - xs0) < 1e-12) xs1 = xs0 + 1.0;

        for (int i = 0; i <= nsamp; ++i) {
            const double t = (double)i / (double)nsamp;
            const double xv = xs0 + t * (xs1 - xs0);
            const double fv = legval(xv);
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
inline double HistoGUIadRoot::legval(double theta)
{
    // Same normalized form as your X11 version
    const double norm = (std::fabs(A0) > 1e-14) ? A0 : 1.0;

    const double c  = std::cos(theta);
    const double c2 = c * c;
    const double c4 = c2 * c2;

    const double aaa = A0 / norm;
    const double aab = (A2E / norm) * (1.5 * c2 - 0.5);
    const double aac = (A4E / norm) * (35.0/8.0 * c4 - 30.0/8.0 * c2 + 3.0/8.0);

    return (aaa + aab + aac);
}

//=============================
// Fit drawing (red only)
//=============================
inline int HistoGUIadRoot::Draw_Fit(double x_low_win, double /*y_low_win*/,
                                    double x_hi_win,  double /*y_hi_win*/)
{
    if (!fit_exists || !fECanvas || !fECanvas->GetCanvas()) return 1;

    // Build fit graph over visible x-range
    const double xl = x_low_win;
    const double xh = x_hi_win;
    if (!(xh > xl)) return 1;

    delete fFitGraph;
    fFitGraph = nullptr;

    const int n = 600;
    fFitGraph = new TGraph(n);

    for (int i = 0; i < n; ++i) {
        const double t = (double)i / (double)(n - 1);
        const double xv = xl + t * (xh - xl);
        const double yv = legval(xv);
        fFitGraph->SetPoint(i, xv, yv);
    }

    fFitGraph->SetLineColor(kRed);
    fFitGraph->SetLineWidth(3);
    fFitGraph->SetMarkerStyle(1);

    // Draw on current pad (assumes frame already drawn)
    fFitGraph->Draw("L SAME");
    return 1;
}

//=============================
// Main data drawing
//=============================
inline int HistoGUIadRoot::DrawData(double x_low_win, double y_low_win,
                                    double x_hi_win,  double y_hi_win)
{
    if (!fECanvas || !fECanvas->GetCanvas()) return 0;

    TCanvas* c = fECanvas->GetCanvas();
    c->cd();
    c->Clear();
    ApplyPadMargins();

    if (!HasUsableData()) {
        // Draw empty frame
        fFrameHist = c->DrawFrame(0.0, 0.0, TMath::Pi(), 1.0);
        StyleFrameHistogram(fFrameHist);
        c->Modified();
        c->Update();
        return 1;
    }

    BuildOrUpdateGraphs(); // ensures graph is synchronized with x/y/ey

    double xl, yl, xh, yh;
    if (x_low_win == -1.0 && y_low_win == -1.0 && x_hi_win == -1.0 && y_hi_win == -1.0) {
        ComputeDefaultBounds(xl, yl, xh, yh);
    } else {
        xl = x_low_win; yl = y_low_win; xh = x_hi_win; yh = y_hi_win;
    }

    min_x = xl; max_x = xh;
    min_y = yl; max_y = yh;

    old_xl = xl; old_xh = xh;
    old_yl = yl; old_yh = yh;

    // Frame / axes
    fFrameHist = c->DrawFrame(xl, yl, xh, yh);
    StyleFrameHistogram(fFrameHist);

    // Data graph + error bars (black)
    if (fGraph) {
        // If you want original normalization behavior (ey/A0 when A0 is set):
        // create a temp graph copy only when A0 != 0. For now, raw ey plotted.
        fGraph->Draw("P E1 SAME");
    }

    // Fit overlay (red)
    if (fit_exists) {
        Draw_Fit(xl, yl, xh, yh);
    }

    // Legend
    delete fLegend;
    fLegend = new TLegend(0.68, 0.78, 0.93, 0.92);
    fLegend->SetBorderSize(1);
    fLegend->SetFillStyle(0);
    fLegend->SetTextColor(kBlack);
    if (fGraph)    fLegend->AddEntry(fGraph, "Data", "pe");
    if (fit_exists && fFitGraph) fLegend->AddEntry(fFitGraph, "Fit", "l");
    fLegend->Draw();

    // Add small annotation for GUI inputs (j1,j2,Egamma,mode)
    {
        // Use ROOT title string only to avoid extra object ownership complexity here
        std::ostringstream title2;
        title2 << plot_title
               << "   [Mode " << fOptionMode << "]"
               << "   j1=" << GetJ1()
               << "   j2=" << GetJ2()
               << "   E#gamma=" << GetGammaKeV() << " keV";
        if (fFrameHist) fFrameHist->SetTitle(title2.str().c_str());
    }

    c->Modified();
    c->Update();
    UpdateTitleText();

    return 1;
}

//=============================
// ROOT GUI event handling
//=============================
inline Bool_t HistoGUIadRoot::ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t /*parm2*/)
{
    switch (GET_MSG(msg)) {
        case kC_COMMAND:
            switch (GET_SUBMSG(msg)) {
                case kCM_BUTTON:
                    switch ((int)parm1) {
                        case kID_Browse: {
                            static TString dir(".");
                            TGFileInfo fi;
                            const char* filetypes[] = {
                                "Data files", "*.csv *.txt *.dat",
                                "CSV files",  "*.csv",
                                "Text files", "*.txt",
                                "DAT files",  "*.dat",
                                "All files",  "*",
                                nullptr, nullptr
                            };
                            fi.fFileTypes = filetypes;
                            fi.fIniDir = StrDup(dir.Data());

                            new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);

                            if (fi.fFilename) {
                                SetFileName(fi.fFilename);
                                dir = fi.fIniDir ? fi.fIniDir : ".";
                            }
                            break;
                        }

                        case kID_Load: {
                            LoadFromGuiFileAndDraw();
                            break;
                        }

                        case kID_Reset: {
                            ResetView();
                            break;
                        }

                        case kID_Quit: {
                            Close();
                            break;
                        }

                        case kID_Opt1: {
                            fOptionMode = 1;
                            UpdateModeButtonStates();
                            UpdateTitleText();
                            // Placeholder for your physics-specific behavior
                            // e.g., recompute coefficients or constraints using j1/j2/Egamma
                            break;
                        }

                        case kID_Opt2: {
                            fOptionMode = 2;
                            UpdateModeButtonStates();
                            UpdateTitleText();
                            // Placeholder for your physics-specific behavior
                            break;
                        }

                        case kID_Opt3: {
                            fOptionMode = 3;
                            UpdateModeButtonStates();
                            UpdateTitleText();
                            // Placeholder for your physics-specific behavior
                            break;
                        }

                        case kID_Redraw: {
                            // Redraw using current j1/j2/Egamma + current data/fit
                            DrawData(-1, -1, -1, -1);
                            break;
                        }

                        default:
                            break;
                    }
                    break;

                default:
                    break;
            }
            break;

        default:
            break;
    }

    return kTRUE;
}

#endif // GUI_AD_ROOT_H
