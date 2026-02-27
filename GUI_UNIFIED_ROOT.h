#ifndef GUI_UNIFIED_ROOT_H
#define GUI_UNIFIED_ROOT_H

//#############################################################################
//
// Unified ROOT GUI (single window) for:
//   (1) Angular Distribution + Fit
//   (2) Chi-squared Scan
//   (3) Close GUI
//
// Improvements:
//   - larger GUI fonts
//   - left-click crosshair with terminal x/y printout
//   - angular plot defaults to [0, pi]
//   - fit always drawn across [0, pi]
//   - chi2 minimum reported automatically
//   - run summary/report box updated every compute
//
//#############################################################################

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>

// ROOT core/graphics
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TString.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVirtualPad.h>
#include <Buttons.h>

// ROOT GUI
#include <RQ_OBJECT.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGLabel.h>
#include <TGLayout.h>
#include <TGTextEntry.h>
#include <TGTextView.h>
#include <TGNumberEntry.h>
#include <TGFileDialog.h>
#include <TRootEmbeddedCanvas.h>
#include <TGResourcePool.h>

class HistoGUIUnifiedRoot : public TGMainFrame {
	RQ_OBJECT("HistoGUIUnifiedRoot")

	public:
		enum EViewMode {
			kViewAngular = 1,
			kViewChi2    = 2
		};

		struct RunRequest {
			double j1 = 0.0;
			double j2 = 0.0;
			double gamma_keV = 0.0;
			double sigma = 0.0;
			std::string dataFile;
		};

		struct PlotResults {
			// Angular
			std::vector<double> dangler;   // degs
			std::vector<double> dydatas;   // normalized y (y/A0)
			std::vector<double> deydata;   // raw yerr
			double fitA0 = 1.0;
			double fitA2 = 0.0;
			double fitA4 = 0.0;
			bool hasFit = false;

			// Chi2 scan
			std::vector<double> tdelta;    // atan(delta)
			std::vector<double> chisqr;    // log(chi^2)

			// Echo inputs
			double j1 = 0.0;
			double j2 = 0.0;
			double gamma_keV = 0.0;
			double sigma = 0.0;

			// Detector attenuation coeffs
			double qd2 = 0.0;
			double qd4 = 0.0;
			bool hasQD = true;
		};

	public:
		HistoGUIUnifiedRoot(const TGWindow* p = gClient->GetRoot(),
				UInt_t w = 1400,
				UInt_t h = 980);

		virtual ~HistoGUIUnifiedRoot();

		int  Init();
		void Close();
		void CloseWindow() override;

		void SetDetectorGeometry(double radius, double distance, double thickness);

		double GetJ1() const;
		double GetJ2() const;
		double GetGammaKeV() const;
		double GetSigma() const;

		void SetJ1(double v);
		void SetJ2(double v);
		void SetGammaKeV(double v);
		void SetSigma(double v);

		void SetDataFileName(const std::string& s);
		std::string GetDataFileName() const;

		void SetViewMode(int mode);
		int  GetViewMode() const { return fViewMode; }

		bool ConsumeRunRequest(RunRequest& out);
		bool ShouldQuitProgram() const { return fQuitRequested; }

		void LoadResults(const PlotResults& r);

		void SetStatus(const std::string& msg);

		void DrawCurrentView(bool autoscale = true);
		void ResetView();
		void Redraw();

		Bool_t ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t parm2) override;

		double legval(double theta) const;
		void PollCanvasClick();
		// void HandleCanvasEvent(Int_t event, Int_t px, Int_t py, TObject* selected);
		void HandleCanvasEvent(int event, int px, int py, TObject* selected);
	private:
		enum EWidgetIds {
			kID_BrowseData   = 1001,
			kID_RunCompute   = 1002,

			kID_Mode1        = 1101,
			kID_Mode2        = 1102,
			kID_Mode3Close   = 1103,

			kID_Redraw       = 1110,
			kID_Reset        = 1111
		};

		void BuildGUI();
		void UpdateWindowTitle();
		void UpdateModeButtonStates();
		void ApplyPadMargins();
		void DrawPlaceholder(const char* msg);

		void DrawAngular(bool autoscale);
		void DrawChi2(bool autoscale);

		bool HasAngularData() const;
		bool HasChi2Data() const;

		void ComputeAngularBounds(double& xl, double& yl, double& xh, double& yh) const;
		void ComputeChi2Bounds(double& xl, double& yl, double& xh, double& yh) const;

		void RebuildAngularGraph();
		void RebuildAngularFitGraph(double xl, double xh);
		void RebuildChi2Graph();

		void StyleFrame(TH1* h, const std::string& title,
				const std::string& xlab,
				const std::string& ylab) const;

		void ApplyUIFont();
		void UpdateReportBox(const PlotResults& r);
		void PrintRunSummaryToTerminal(const PlotResults& r) const;

		void ClearCrosshair();
		void DrawCrosshairOverlay();
		bool GetCurrentPlotUserCoordsFromPixel(Int_t px, Int_t py, double& x, double& y) const;

	private:
		//================ GUI widgets ================
		TGVerticalFrame*     fMainV;

		TGHorizontalFrame*   fRowGeom;
		TGHorizontalFrame*   fRowParams;
		TGHorizontalFrame*   fRowFile;
		TGHorizontalFrame*   fRowButtons;
		TGHorizontalFrame*   fRowReport;
		TGHorizontalFrame*   fRowStatus;

		TGLabel*             fLblGeom;
		TGLabel*             fLblGeomVal;

		TGLabel*             fLblJ1;
		TGLabel*             fLblJ2;
		TGLabel*             fLblEg;
		TGLabel*             fLblSigma;

		TGNumberEntry*       fEntJ1;
		TGNumberEntry*       fEntJ2;
		TGNumberEntry*       fEntEg;
		TGNumberEntry*       fEntSigma;

		TGLabel*             fLblDataFile;
		TGTextEntry*         fEntDataFile;
		TGTextButton*        fBtnBrowseData;
		TGTextButton*        fBtnRun;

		TGTextButton*        fBtnMode1;
		TGTextButton*        fBtnMode2;
		TGTextButton*        fBtnMode3;
		TGTextButton*        fBtnRedraw;
		TGTextButton*        fBtnReset;

		TGTextView*          fReportView;
		TGLabel*             fLblStatus;

		TRootEmbeddedCanvas* fECanvas;

		// GUI font
		FontStruct_t         fUIFont;
		bool                 fCanvasConnected;

		//================ State ================
		int    fViewMode;
		bool   fIsMapped;
		bool   fQuitRequested;

		bool       fRunRequested;
		RunRequest fPendingReq;

		bool   fAngularHasLastRange;
		bool   fChi2HasLastRange;

		double fAng_xl, fAng_xh, fAng_yl, fAng_yh;
		double fChi_xl, fChi_xh, fChi_yl, fChi_yh;

		double fMarginLeft;
		double fMarginRight;
		double fMarginTop;
		double fMarginBottom;

		// detector geometry
		double fDetRadius;
		double fDetDistance;
		double fDetThickness;

		//================ Data / plot objects ================
		std::vector<double> fAngX;
		std::vector<double> fAngY;
		std::vector<double> fAngEY;

		bool   fFitExists;
		double A0;
		double A2E;
		double A4E;

		std::vector<double> fChi2X;
		std::vector<double> fChi2Y;

		// last run summary values
		double fLastBestAtanDelta;
		double fLastBestDelta;
		double fLastBestLogChi2;
		double fLastBestChi2;
		bool   fHasBestChi2;

		double fLastQD2;
		double fLastQD4;
		bool   fHasQD;

		std::string fAngularTitle;
		std::string fAngularXLabel;
		std::string fAngularYLabel;

		std::string fChi2Title;
		std::string fChi2XLabel;
		std::string fChi2YLabel;

		TGraphErrors* fGraphAngular;
		TGraph*       fGraphAngFit;
		TGraph*       fGraphChi2;
		TLegend*      fLegend;
		TH1*          fFrameHist; // owned by pad

		// Crosshair overlay

		bool          fCrosshairActive;
		int           fCrosshairViewMode;
		double        fCrossX;
		double        fCrossY;
		TLine*        fCrossVLine;
		TLine*        fCrossHLine;

		int           fLastCanvasEvent;
		int           fLastCanvasPX;
		int           fLastCanvasPY;

};

//==============================================================================
// Constructor / Destructor
//==============================================================================
inline HistoGUIUnifiedRoot::HistoGUIUnifiedRoot(const TGWindow* p, UInt_t w, UInt_t h)
	: TGMainFrame(p, w, h),
	fMainV(nullptr),
	fRowGeom(nullptr),
	fRowParams(nullptr),
	fRowFile(nullptr),
	fRowButtons(nullptr),
	fRowReport(nullptr),
	fRowStatus(nullptr),
	fLblGeom(nullptr),
	fLblGeomVal(nullptr),
	fLblJ1(nullptr), fLblJ2(nullptr), fLblEg(nullptr), fLblSigma(nullptr),
	fEntJ1(nullptr), fEntJ2(nullptr), fEntEg(nullptr), fEntSigma(nullptr),
	fLblDataFile(nullptr), fEntDataFile(nullptr), fBtnBrowseData(nullptr), fBtnRun(nullptr),
	fBtnMode1(nullptr), fBtnMode2(nullptr), fBtnMode3(nullptr),
	fBtnRedraw(nullptr), fBtnReset(nullptr),
	fReportView(nullptr),
	fLblStatus(nullptr),
	fECanvas(nullptr),
	fUIFont(0),
	fCanvasConnected(false),
	fViewMode(kViewAngular),
	fIsMapped(false),
	fQuitRequested(false),
	fRunRequested(false),
	fAngularHasLastRange(false),
	fChi2HasLastRange(false),
	fAng_xl(0.0), fAng_xh(1.0), fAng_yl(0.0), fAng_yh(1.0),
	fChi_xl(0.0), fChi_xh(1.0), fChi_yl(0.0), fChi_yh(1.0),
	fMarginLeft(0.12), fMarginRight(0.04), fMarginTop(0.08), fMarginBottom(0.12),
	fDetRadius(0.0), fDetDistance(0.0), fDetThickness(0.0),
	fFitExists(false),
	A0(1.0), A2E(0.0), A4E(0.0),
	fLastBestAtanDelta(0.0), fLastBestDelta(0.0), fLastBestLogChi2(0.0), fLastBestChi2(0.0),
	fHasBestChi2(false),
	fLastQD2(0.0), fLastQD4(0.0), fHasQD(false),
	fAngularTitle("Angular Distribution + Fit"),
	fAngularXLabel("Theta (deg)"),
	fAngularYLabel("W(#theta) / A_{0}"),
	fChi2Title("Chi-Squared Scan"),
	fChi2XLabel("atan(#delta) [rad]"),
	fChi2YLabel("log(#chi^{2})"),
	fGraphAngular(nullptr),
	fGraphAngFit(nullptr),
	fGraphChi2(nullptr),
	fLegend(nullptr),
	fFrameHist(nullptr),
	fCrosshairActive(false),
	fCrosshairViewMode(-1),
	fCrossX(0.0),
	fCrossY(0.0),
	fCrossVLine(nullptr),
	fCrossHLine(nullptr),
	fLastCanvasEvent(-9999),
	fLastCanvasPX(-9999),
	fLastCanvasPY(-9999)
{
	SetCleanup(kDeepCleanup);
	BuildGUI();
	ApplyUIFont();
	UpdateModeButtonStates();
	UpdateWindowTitle();
	SetStatus("Ready.");
}

inline HistoGUIUnifiedRoot::~HistoGUIUnifiedRoot()
{
	delete fGraphAngular; fGraphAngular = nullptr;
	delete fGraphAngFit;  fGraphAngFit  = nullptr;
	delete fGraphChi2;    fGraphChi2    = nullptr;
	delete fLegend;       fLegend       = nullptr;
	delete fCrossVLine;   fCrossVLine   = nullptr;
	delete fCrossHLine;   fCrossHLine   = nullptr;
	Cleanup();
}

//==============================================================================
// GUI builders
//==============================================================================
inline void HistoGUIUnifiedRoot::BuildGUI()
{
	fMainV = new TGVerticalFrame(this);

	// ---------------- Row: detector geometry display ----------------
	fRowGeom = new TGHorizontalFrame(fMainV);
	fLblGeom = new TGLabel(fRowGeom, "Detector geometry (from AD3.cxx):");
	fLblGeomVal = new TGLabel(fRowGeom, "R=?, D=?, T=?");
	fRowGeom->AddFrame(fLblGeom,    new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 6, 4, 4));
	fRowGeom->AddFrame(fLblGeomVal, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 6, 4, 4));
	fRowGeom->AddFrame(new TGLabel(fRowGeom, " "), new TGLayoutHints(kLHintsExpandX));

	// ---------------- Row: physics params ----------------
	fRowParams = new TGHorizontalFrame(fMainV);

	fLblJ1 = new TGLabel(fRowParams, "j1:");
	fEntJ1 = new TGNumberEntry(fRowParams, 0.0, 8, -1, TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber);

	fLblJ2 = new TGLabel(fRowParams, "j2:");
	fEntJ2 = new TGNumberEntry(fRowParams, 0.0, 8, -1, TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber);

	fLblEg = new TGLabel(fRowParams, "E#gamma (keV):");
	fEntEg = new TGNumberEntry(fRowParams, 0.0, 10, -1, TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber);

	fLblSigma = new TGLabel(fRowParams, "Sigma:");
	fEntSigma = new TGNumberEntry(fRowParams, 0.0, 8, -1, TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber);

	fRowParams->AddFrame(fLblJ1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 4, 4));
	fRowParams->AddFrame(fEntJ1, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 10, 4, 4));

	fRowParams->AddFrame(fLblJ2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 4, 4));
	fRowParams->AddFrame(fEntJ2, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 10, 4, 4));

	fRowParams->AddFrame(fLblEg, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 4, 4));
	fRowParams->AddFrame(fEntEg, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 10, 4, 4));

	fRowParams->AddFrame(fLblSigma, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 4, 4));
	fRowParams->AddFrame(fEntSigma, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 6, 4, 4));

	fRowParams->AddFrame(new TGLabel(fRowParams, " "), new TGLayoutHints(kLHintsExpandX));

	// ---------------- Row: data file + Run ----------------
	fRowFile = new TGHorizontalFrame(fMainV);
	fLblDataFile   = new TGLabel(fRowFile, "Angular data file:");
	fEntDataFile   = new TGTextEntry(fRowFile, "");
	fBtnBrowseData = new TGTextButton(fRowFile, "Browse", kID_BrowseData);
	fBtnRun        = new TGTextButton(fRowFile, "Run / Compute", kID_RunCompute);

	fBtnBrowseData->Associate(this);
	fBtnRun->Associate(this);

	fRowFile->AddFrame(fLblDataFile,   new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 4, 4));
	fRowFile->AddFrame(fEntDataFile,   new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 6, 4, 4));
	fRowFile->AddFrame(fBtnBrowseData, new TGLayoutHints(kLHintsRight | kLHintsCenterY, 2, 4, 4, 4));
	fRowFile->AddFrame(fBtnRun,        new TGLayoutHints(kLHintsRight | kLHintsCenterY, 2, 6, 4, 4));

	// ---------------- Row: mode buttons + actions ----------------
	fRowButtons = new TGHorizontalFrame(fMainV);

	fBtnMode1  = new TGTextButton(fRowButtons, "1", kID_Mode1);
	fBtnMode2  = new TGTextButton(fRowButtons, "2", kID_Mode2);
	fBtnMode3  = new TGTextButton(fRowButtons, "3", kID_Mode3Close);

	fBtnRedraw = new TGTextButton(fRowButtons, "Redraw", kID_Redraw);
	fBtnReset  = new TGTextButton(fRowButtons, "Reset Zoom", kID_Reset);

	fBtnMode1->Associate(this);
	fBtnMode2->Associate(this);
	fBtnMode3->Associate(this);
	fBtnRedraw->Associate(this);
	fBtnReset->Associate(this);

	fRowButtons->AddFrame(fBtnMode1,  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 2, 4, 4));
	fRowButtons->AddFrame(fBtnMode2,  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 4, 4));
	fRowButtons->AddFrame(fBtnMode3,  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 12, 4, 4));
	fRowButtons->AddFrame(fBtnRedraw, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 4, 4, 4));
	fRowButtons->AddFrame(fBtnReset,  new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 4, 4, 4));
	fRowButtons->AddFrame(new TGLabel(fRowButtons, " "), new TGLayoutHints(kLHintsExpandX));

	// ---------------- Row: report summary ----------------
	fRowReport = new TGHorizontalFrame(fMainV);
	fReportView = new TGTextView(fRowReport, 1200, 120);
	fReportView->LoadBuffer(
			"Run summary:\n"
			"  No results yet.\n"
			);
	fRowReport->AddFrame(fReportView, new TGLayoutHints(kLHintsExpandX, 6, 6, 2, 4));

	// Embedded canvas
	fECanvas = new TRootEmbeddedCanvas("unified_canvas", fMainV, 1000, 700);

	// Status row
	fRowStatus = new TGHorizontalFrame(fMainV);
	fLblStatus = new TGLabel(fRowStatus, "Status: ");
	fRowStatus->AddFrame(fLblStatus, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 6, 4, 4));
	fRowStatus->AddFrame(new TGLabel(fRowStatus, " "), new TGLayoutHints(kLHintsExpandX));

	// Assemble
	fMainV->AddFrame(fRowGeom,    new TGLayoutHints(kLHintsTop | kLHintsExpandX));
	fMainV->AddFrame(fRowParams,  new TGLayoutHints(kLHintsTop | kLHintsExpandX));
	fMainV->AddFrame(fRowFile,    new TGLayoutHints(kLHintsTop | kLHintsExpandX));
	fMainV->AddFrame(fRowButtons, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
	fMainV->AddFrame(fRowReport,  new TGLayoutHints(kLHintsTop | kLHintsExpandX));
	fMainV->AddFrame(fECanvas,    new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 6, 6, 4, 6));
	fMainV->AddFrame(fRowStatus,  new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

	AddFrame(fMainV, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

	SetWindowName("Unified Angular/Chi2 GUI (ROOT)");
	MapSubwindows();
	Resize(GetDefaultSize());

	if (fECanvas && fECanvas->GetCanvas()) {
		TCanvas* c = fECanvas->GetCanvas();
		c->SetFillColor(0);
		c->SetBorderMode(0);
		c->SetFrameBorderMode(0);
		c->SetGrid(0, 0);
	}
}


inline void HistoGUIUnifiedRoot::ApplyUIFont()
{
	// Larger GUI font that works with this ROOT API.
	fUIFont = gClient->GetFontByName("-adobe-helvetica-medium-r-*-*-16-*-*-*-*-*-iso8859-1");
	if (!fUIFont) {
		fUIFont = gClient->GetFontByName("-*-helvetica-medium-r-*-*-16-*-*-*-*-*-iso8859-1");
	}
	if (!fUIFont) return;

	// Labels
	if (fLblGeom)      fLblGeom->SetTextFont(fUIFont, kFALSE);
	if (fLblGeomVal)   fLblGeomVal->SetTextFont(fUIFont, kFALSE);
	if (fLblJ1)        fLblJ1->SetTextFont(fUIFont, kFALSE);
	if (fLblJ2)        fLblJ2->SetTextFont(fUIFont, kFALSE);
	if (fLblEg)        fLblEg->SetTextFont(fUIFont, kFALSE);
	if (fLblSigma)     fLblSigma->SetTextFont(fUIFont, kFALSE);
	if (fLblDataFile)  fLblDataFile->SetTextFont(fUIFont, kFALSE);
	if (fLblStatus)    fLblStatus->SetTextFont(fUIFont, kFALSE);

	// Buttons
	if (fBtnBrowseData) fBtnBrowseData->SetFont(fUIFont, kFALSE);
	if (fBtnRun)        fBtnRun->SetFont(fUIFont, kFALSE);
	if (fBtnMode1)      fBtnMode1->SetFont(fUIFont, kFALSE);
	if (fBtnMode2)      fBtnMode2->SetFont(fUIFont, kFALSE);
	if (fBtnMode3)      fBtnMode3->SetFont(fUIFont, kFALSE);
	if (fBtnRedraw)     fBtnRedraw->SetFont(fUIFont, kFALSE);
	if (fBtnReset)      fBtnReset->SetFont(fUIFont, kFALSE);

	// Text entry
	if (fEntDataFile)   fEntDataFile->SetFont(fUIFont, kTRUE);

	// Report box
	if (fReportView)    fReportView->SetFont(fUIFont);

	// Number-entry text fields
	if (fEntJ1 && fEntJ1->GetNumberEntry())         fEntJ1->GetNumberEntry()->SetFont(fUIFont, kTRUE);
	if (fEntJ2 && fEntJ2->GetNumberEntry())         fEntJ2->GetNumberEntry()->SetFont(fUIFont, kTRUE);
	if (fEntEg && fEntEg->GetNumberEntry())         fEntEg->GetNumberEntry()->SetFont(fUIFont, kTRUE);
	if (fEntSigma && fEntSigma->GetNumberEntry())   fEntSigma->GetNumberEntry()->SetFont(fUIFont, kTRUE);

	// Make widgets taller/wider so the larger font actually fits visually.
	if (fEntJ1)      fEntJ1->Resize(110, 30);
	if (fEntJ2)      fEntJ2->Resize(110, 30);
	if (fEntEg)      fEntEg->Resize(130, 30);
	if (fEntSigma)   fEntSigma->Resize(110, 30);

	if (fEntDataFile)    fEntDataFile->Resize(520, 30);

	if (fBtnBrowseData)  fBtnBrowseData->Resize(110, 32);
	if (fBtnRun)         fBtnRun->Resize(170, 32);

	if (fBtnMode1)       fBtnMode1->Resize(50, 32);
	if (fBtnMode2)       fBtnMode2->Resize(50, 32);
	if (fBtnMode3)       fBtnMode3->Resize(50, 32);
	if (fBtnRedraw)      fBtnRedraw->Resize(110, 32);
	if (fBtnReset)       fBtnReset->Resize(130, 32);

	if (fReportView)     fReportView->Resize(1200, 135);

	Layout();
}


inline void HistoGUIUnifiedRoot::ApplyPadMargins()
{
	if (!fECanvas || !fECanvas->GetCanvas()) return;
	TCanvas* c = fECanvas->GetCanvas();
	c->SetLeftMargin(fMarginLeft);
	c->SetRightMargin(fMarginRight);
	c->SetTopMargin(fMarginTop);
	c->SetBottomMargin(fMarginBottom);
}

inline void HistoGUIUnifiedRoot::UpdateModeButtonStates()
{
	if (fBtnMode1) fBtnMode1->SetText(fViewMode == kViewAngular ? "[1]" : "1");
	if (fBtnMode2) fBtnMode2->SetText(fViewMode == kViewChi2    ? "[2]" : "2");
	if (fBtnMode3) fBtnMode3->SetText("3");
}

inline void HistoGUIUnifiedRoot::UpdateWindowTitle()
{
	std::ostringstream os;
	os << "Unified Angular/Chi2 GUI (ROOT)"
		<< " | View=" << (fViewMode == kViewAngular ? "Angular" : "Chi2")
		<< " | j1=" << GetJ1()
		<< " | j2=" << GetJ2()
		<< " | E#gamma=" << GetGammaKeV() << " keV"
		<< " | sigma=" << GetSigma();
	SetWindowName(os.str().c_str());
}

inline void HistoGUIUnifiedRoot::StyleFrame(TH1* h,
		const std::string& title,
		const std::string& xlab,
		const std::string& ylab) const
{
	if (!h) return;

	h->SetTitle(title.c_str());
	h->SetStats(0);
	h->SetLineColor(kBlack);

	h->GetXaxis()->SetTitle(xlab.c_str());
	h->GetYaxis()->SetTitle(ylab.c_str());

	h->GetXaxis()->SetTitleSize(0.055);
	h->GetYaxis()->SetTitleSize(0.055);
	h->GetXaxis()->SetLabelSize(0.045);
	h->GetYaxis()->SetLabelSize(0.045);
	h->GetXaxis()->SetTitleOffset(1.00);
	h->GetYaxis()->SetTitleOffset(1.00);
}
inline void HistoGUIUnifiedRoot::PollCanvasClick()
{
	if (!fECanvas || !fECanvas->GetCanvas()) return;
	if (!fFrameHist) return;

	TCanvas* c = fECanvas->GetCanvas();

	const int event = c->GetEvent();
	const int px    = c->GetEventX();
	const int py    = c->GetEventY();

	// Ignore unchanged event state so we do not retrigger forever
	if (event == fLastCanvasEvent &&
			px    == fLastCanvasPX &&
			py    == fLastCanvasPY) {
		return;
	}

	fLastCanvasEvent = event;
	fLastCanvasPX    = px;
	fLastCanvasPY    = py;

	if (event != kButton1Down) return;

	HandleCanvasEvent(event, px, py, nullptr);
}
//==============================================================================
// Lifecycle
//==============================================================================
inline int HistoGUIUnifiedRoot::Init()
{
	if (!fIsMapped) {
		MapWindow();
		fIsMapped = true;
	}

	DrawCurrentView(true);
	return 1;
}
inline void HistoGUIUnifiedRoot::Close()
{
	fQuitRequested = true;
	UnmapWindow();
}

inline void HistoGUIUnifiedRoot::CloseWindow()
{
	Close();
}

//==============================================================================
// Detector geometry
//==============================================================================
inline void HistoGUIUnifiedRoot::SetDetectorGeometry(double radius, double distance, double thickness)
{
	fDetRadius = radius;
	fDetDistance = distance;
	fDetThickness = thickness;

	if (fLblGeomVal) {
		std::ostringstream os;
		os << "R=" << fDetRadius << "  D=" << fDetDistance << "  T=" << fDetThickness;
		fLblGeomVal->SetText(os.str().c_str());
		Layout();
	}
}

//==============================================================================
// Shared GUI parameters
//==============================================================================
inline double HistoGUIUnifiedRoot::GetJ1() const        { return fEntJ1 ? fEntJ1->GetNumber() : 0.0; }
inline double HistoGUIUnifiedRoot::GetJ2() const        { return fEntJ2 ? fEntJ2->GetNumber() : 0.0; }
inline double HistoGUIUnifiedRoot::GetGammaKeV() const  { return fEntEg ? fEntEg->GetNumber() : 0.0; }
inline double HistoGUIUnifiedRoot::GetSigma() const     { return fEntSigma ? fEntSigma->GetNumber() : 0.0; }

inline void HistoGUIUnifiedRoot::SetJ1(double v)        { if (fEntJ1) fEntJ1->SetNumber(v); UpdateWindowTitle(); }
inline void HistoGUIUnifiedRoot::SetJ2(double v)        { if (fEntJ2) fEntJ2->SetNumber(v); UpdateWindowTitle(); }
inline void HistoGUIUnifiedRoot::SetGammaKeV(double v)  { if (fEntEg) fEntEg->SetNumber(v); UpdateWindowTitle(); }
inline void HistoGUIUnifiedRoot::SetSigma(double v)     { if (fEntSigma) fEntSigma->SetNumber(v); UpdateWindowTitle(); }

//==============================================================================
// File path
//==============================================================================
inline void HistoGUIUnifiedRoot::SetDataFileName(const std::string& s)
{
	if (fEntDataFile) fEntDataFile->SetText(s.c_str());
}

inline std::string HistoGUIUnifiedRoot::GetDataFileName() const
{
	return (fEntDataFile && fEntDataFile->GetText()) ? std::string(fEntDataFile->GetText()) : std::string();
}

//==============================================================================
// View control
//==============================================================================
inline void HistoGUIUnifiedRoot::SetViewMode(int mode)
{
	if (mode != kViewAngular && mode != kViewChi2) return;
	fViewMode = mode;
	UpdateModeButtonStates();
	UpdateWindowTitle();
	DrawCurrentView(true);
}

//==============================================================================
// Run/Compute handshake
//==============================================================================
inline bool HistoGUIUnifiedRoot::ConsumeRunRequest(RunRequest& out)
{
	if (!fRunRequested) return false;
	out = fPendingReq;
	fRunRequested = false;
	return true;
}

//==============================================================================
// Status
//==============================================================================
inline void HistoGUIUnifiedRoot::SetStatus(const std::string& msg)
{
	if (!fLblStatus) return;
	std::string s = "Status: " + msg;
	fLblStatus->SetText(s.c_str());
	Layout();
}

//==============================================================================
// Results ingestion
//==============================================================================
inline void HistoGUIUnifiedRoot::LoadResults(const PlotResults& r)
{
	std::printf("LoadResults called.\n");

	fAngX  = r.dangler;
	fAngY  = r.dydatas;
	fAngEY = r.deydata;

	fChi2X = r.tdelta;
	fChi2Y = r.chisqr;

	fFitExists = r.hasFit;
	if (fFitExists) {
		A0  = r.fitA0;
		A2E = r.fitA2;
		A4E = r.fitA4;
	} else {
		A0 = 1.0; A2E = 0.0; A4E = 0.0;
	}

	//  fHasQD = r.hasQD;
	// fLastQD2 = r.qd2;
	// fLastQD4 = r.qd4;

	fLastQD2 = r.qd2;
	fLastQD4 = r.qd4;
	fHasQD   = std::isfinite(fLastQD2) && std::isfinite(fLastQD4);
	// Chi2 minimum
	fHasBestChi2 = false;
	fLastBestAtanDelta = 0.0;
	fLastBestDelta = 0.0;
	fLastBestLogChi2 = 0.0;
	fLastBestChi2 = 0.0;

	if (!fChi2X.empty() && fChi2X.size() == fChi2Y.size()) {
		size_t imin = 0;
		for (size_t i = 1; i < fChi2Y.size(); ++i) {
			if (fChi2Y[i] < fChi2Y[imin]) imin = i;
		}
		fLastBestAtanDelta = fChi2X[imin];
		fLastBestDelta     = std::tan(fLastBestAtanDelta);
		fLastBestLogChi2   = fChi2Y[imin];
		fLastBestChi2      = std::exp(fLastBestLogChi2);
		fHasBestChi2       = true;
	}

	RebuildAngularGraph();
	RebuildChi2Graph();

	UpdateReportBox(r);
	PrintRunSummaryToTerminal(r);

	ClearCrosshair();
	DrawCurrentView(true);
}

//==============================================================================
// Drawing helpers
//==============================================================================
inline bool HistoGUIUnifiedRoot::HasAngularData() const { return (!fAngX.empty() && fAngX.size() == fAngY.size()); }
inline bool HistoGUIUnifiedRoot::HasChi2Data() const    { return (!fChi2X.empty() && fChi2X.size() == fChi2Y.size()); }

inline void HistoGUIUnifiedRoot::RebuildAngularGraph()
{
	delete fGraphAngular;
	fGraphAngular = nullptr;

	if (!HasAngularData()) return;

	const int n = (int)fAngX.size();
	fGraphAngular = new TGraphErrors(n);

	for (int i = 0; i < n; ++i) {
		double ey = (i < (int)fAngEY.size()) ? fAngEY[(size_t)i] : 0.0;
		if (fFitExists && std::fabs(A0) > 1e-12) ey /= A0;

		const double xdeg = fAngX[(size_t)i] * 180.0 / TMath::Pi();

		fGraphAngular->SetPoint(i, xdeg, fAngY[(size_t)i]);
		fGraphAngular->SetPointError(i, 0.0, ey);
	}

	fGraphAngular->SetLineColor(kBlack);
	fGraphAngular->SetMarkerColor(kBlack);
	fGraphAngular->SetMarkerStyle(20);
	fGraphAngular->SetMarkerSize(1.0);
	fGraphAngular->SetLineWidth(2);
}


inline void HistoGUIUnifiedRoot::RebuildChi2Graph()
{
	delete fGraphChi2;
	fGraphChi2 = nullptr;

	if (!HasChi2Data()) return;

	const int n = (int)fChi2X.size();
	fGraphChi2 = new TGraph(n);

	for (int i = 0; i < n; ++i) {
		fGraphChi2->SetPoint(i, fChi2X[(size_t)i], fChi2Y[(size_t)i]);
	}

	fGraphChi2->SetLineColor(kBlack);
	fGraphChi2->SetMarkerColor(kBlack);
	fGraphChi2->SetMarkerStyle(20);
	fGraphChi2->SetMarkerSize(0.9);
	fGraphChi2->SetLineWidth(2);
}
inline void HistoGUIUnifiedRoot::RebuildAngularFitGraph(double xl_deg, double xh_deg)
{
	delete fGraphAngFit;
	fGraphAngFit = nullptr;

	if (!fFitExists) return;
	if (!(xh_deg > xl_deg)) return;

	const int n = 600;
	fGraphAngFit = new TGraph(n);

	for (int i = 0; i < n; ++i) {
		const double t = (double)i / (double)(n - 1);
		const double xdeg = xl_deg + t * (xh_deg - xl_deg);
		const double xrad = xdeg * TMath::Pi() / 180.0;
		const double yv   = legval(xrad);
		fGraphAngFit->SetPoint(i, xdeg, yv);
	}

	fGraphAngFit->SetLineColor(kRed);
	fGraphAngFit->SetLineWidth(3);
	fGraphAngFit->SetMarkerStyle(1);
}


inline void HistoGUIUnifiedRoot::ComputeAngularBounds(double& xl, double& yl, double& xh, double& yh) const
{
	// Always show angular distribution over 0 -> 180 degrees by default
	xl = 0.0;
	xh = 180.0;

	if (!HasAngularData()) {
		yl = 0.0;
		yh = 1.0;
		return;
	}

	double y_min = fAngY[0];
	double y_max = fAngY[0];

	for (size_t i = 0; i < fAngX.size(); ++i) {
		double ey = (i < fAngEY.size()) ? fAngEY[i] : 0.0;
		if (fFitExists && std::fabs(A0) > 1e-12) ey /= A0;

		y_min = std::min(y_min, fAngY[i] - ey);
		y_max = std::max(y_max, fAngY[i] + ey);
	}

	// Sample fit over full 0 -> 180 degree range
	if (fFitExists) {
		const int nsamp = 600;
		for (int i = 0; i <= nsamp; ++i) {
			const double t    = (double)i / (double)nsamp;
			const double xdeg = xl + t * (xh - xl);
			const double xrad = xdeg * TMath::Pi() / 180.0;
			const double fv   = legval(xrad);
			y_min = std::min(y_min, fv);
			y_max = std::max(y_max, fv);
		}
	}

	double y_pad = 0.10 * (y_max - y_min);
	if (std::fabs(y_pad) < 1e-12) y_pad = 0.5;

	yl = y_min - y_pad;
	yh = y_max + y_pad;
}


inline void HistoGUIUnifiedRoot::ComputeChi2Bounds(double& xl, double& yl, double& xh, double& yh) const
{
	if (!HasChi2Data()) {
		xl = -TMath::Pi()/2.0; xh = TMath::Pi()/2.0;
		yl = 0.0; yh = 1.0;
		return;
	}

	double x_min = fChi2X[0], x_max = fChi2X[0];
	double y_min = fChi2Y[0], y_max = fChi2Y[0];

	for (size_t i = 0; i < fChi2X.size(); ++i) {
		x_min = std::min(x_min, fChi2X[i]);
		x_max = std::max(x_max, fChi2X[i]);
		y_min = std::min(y_min, fChi2Y[i]);
		y_max = std::max(y_max, fChi2Y[i]);
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

inline void HistoGUIUnifiedRoot::DrawPlaceholder(const char* msg)
{
	if (!fECanvas || !fECanvas->GetCanvas()) return;
	TCanvas* c = fECanvas->GetCanvas();
	c->cd();
	c->Clear();
	ApplyPadMargins();

	fFrameHist = c->DrawFrame(0.0, 0.0, 1.0, 1.0);
	StyleFrame(fFrameHist, "No Data", "x", "y");

	TLatex latex;
	latex.SetNDC(true);
	latex.SetTextColor(kBlack);
	latex.SetTextSize(0.04);
	latex.DrawLatex(0.18, 0.55, msg);

	c->Modified();
	c->Update();
}

inline void HistoGUIUnifiedRoot::DrawAngular(bool autoscale)
{
	if (!fECanvas || !fECanvas->GetCanvas()) return;

	TCanvas* c = fECanvas->GetCanvas();
	c->cd();
	c->Clear();
	ApplyPadMargins();

	if (!HasAngularData()) {
		DrawPlaceholder("No angular-distribution results yet. Click Run / Compute.");
		return;
	}

	RebuildAngularGraph();

	double xl, yl, xh, yh;
	if (!autoscale && fAngularHasLastRange) {
		xl = fAng_xl; yl = fAng_yl; xh = fAng_xh; yh = fAng_yh;
	} else {
		ComputeAngularBounds(xl, yl, xh, yh);
	}

	fAng_xl = xl; fAng_yl = yl; fAng_xh = xh; fAng_yh = yh;
	fAngularHasLastRange = true;

	std::ostringstream title;
	title << fAngularTitle
		<< "   [1]"
		<< "   j1=" << GetJ1()
		<< "   j2=" << GetJ2()
		<< "   E#gamma=" << GetGammaKeV() << " keV"
		<< "   #sigma=" << GetSigma();

	fFrameHist = c->DrawFrame(xl, yl, xh, yh);
	StyleFrame(fFrameHist, title.str(), fAngularXLabel, fAngularYLabel);

	if (fGraphAngular) fGraphAngular->Draw("P E1 SAME");

	// Always draw fit over full 0 -> pi range
	if (fFitExists) {
		RebuildAngularFitGraph(0.0, 180.);
		if (fGraphAngFit) fGraphAngFit->Draw("L SAME");
	}

	delete fLegend;
	fLegend = new TLegend(0.68, 0.78, 0.93, 0.92);
	fLegend->SetBorderSize(1);
	fLegend->SetFillStyle(0);
	fLegend->SetTextColor(kBlack);
	if (fGraphAngular) fLegend->AddEntry(fGraphAngular, "Data", "pe");
	if (fFitExists && fGraphAngFit) fLegend->AddEntry(fGraphAngFit, "Fit", "l");
	fLegend->Draw();

	DrawCrosshairOverlay();

	c->Modified();
	c->Update();
}

inline void HistoGUIUnifiedRoot::DrawChi2(bool autoscale)
{
	if (!fECanvas || !fECanvas->GetCanvas()) return;

	TCanvas* c = fECanvas->GetCanvas();
	c->cd();
	c->Clear();
	ApplyPadMargins();

	if (!HasChi2Data()) {
		DrawPlaceholder("No chi-squared results yet. Click Run / Compute.");
		return;
	}

	RebuildChi2Graph();

	double xl, yl, xh, yh;
	if (!autoscale && fChi2HasLastRange) {
		xl = fChi_xl; yl = fChi_yl; xh = fChi_xh; yh = fChi_yh;
	} else {
		ComputeChi2Bounds(xl, yl, xh, yh);
	}

	fChi_xl = xl; fChi_yl = yl; fChi_xh = xh; fChi_yh = yh;
	fChi2HasLastRange = true;

	std::ostringstream title;
	title << fChi2Title
		<< "   [2]"
		<< "   j1=" << GetJ1()
		<< "   j2=" << GetJ2()
		<< "   E#gamma=" << GetGammaKeV() << " keV"
		<< "   #sigma=" << GetSigma();

	fFrameHist = c->DrawFrame(xl, yl, xh, yh);
	StyleFrame(fFrameHist, title.str(), fChi2XLabel, fChi2YLabel);

	if (fGraphChi2) fGraphChi2->Draw("LP SAME");

	delete fLegend;
	fLegend = new TLegend(0.73, 0.83, 0.93, 0.92);
	fLegend->SetBorderSize(1);
	fLegend->SetFillStyle(0);
	fLegend->SetTextColor(kBlack);
	if (fGraphChi2) fLegend->AddEntry(fGraphChi2, "Chi^{2} scan", "lp");
	fLegend->Draw();

	DrawCrosshairOverlay();

	c->Modified();
	c->Update();
}

inline void HistoGUIUnifiedRoot::DrawCurrentView(bool autoscale)
{
	UpdateModeButtonStates();
	UpdateWindowTitle();

	if (fViewMode == kViewAngular) DrawAngular(autoscale);
	else                           DrawChi2(autoscale);
}

inline void HistoGUIUnifiedRoot::ResetView()
{
	if (fViewMode == kViewAngular) fAngularHasLastRange = false;
	if (fViewMode == kViewChi2)    fChi2HasLastRange = false;
	DrawCurrentView(true);
}

inline void HistoGUIUnifiedRoot::Redraw()
{
	DrawCurrentView(true);
}

//==============================================================================
// Angular fit function (normalized)
//==============================================================================
inline double HistoGUIUnifiedRoot::legval(double theta) const
{
	const double norm = (std::fabs(A0) > 1e-14) ? A0 : 1.0;

	const double c  = std::cos(theta);
	const double c2 = c * c;
	const double c4 = c2 * c2;

	const double p2 = (1.5 * c2 - 0.5);
	const double p4 = (35.0/8.0 * c4 - 30.0/8.0 * c2 + 3.0/8.0);

	return (A0 / norm) + (A2E / norm) * p2 + (A4E / norm) * p4;
}

//==============================================================================
// Report helpers
//==============================================================================

inline void HistoGUIUnifiedRoot::UpdateReportBox(const PlotResults& r)
{
	if (!fReportView) return;

	const double j1v = (std::fabs(r.j1) > 1e-14) ? r.j1 : GetJ1();
	const double j2v = (std::fabs(r.j2) > 1e-14) ? r.j2 : GetJ2();
	const double Eg  = (std::fabs(r.gamma_keV) > 1e-14) ? r.gamma_keV : GetGammaKeV();
	const double sig = (std::fabs(r.sigma) > 1e-14) ? r.sigma : GetSigma();

	std::ostringstream os;
	os << std::fixed << std::setprecision(6);
	os << "Run summary\n";
	os << "  Inputs : E_gamma = " << Eg << " keV"
		<< "   j1 = " << j1v
		<< "   j2 = " << j2v
		<< "   sigma = " << sig << "\n";

	os << "  Fit    : A0 = " << A0
		<< "   A2 = " << A2E
		<< "   A4 = " << A4E << "\n";

	if (fHasBestChi2) {
		os << "  Best   : atan(delta) = " << fLastBestAtanDelta << " rad";
	} else {
	}

	if (fHasQD) {
		os << "  QD     : QD2 = " << fLastQD2
			<< "   QD4 = " << fLastQD4 << "\n";
	} else {
		os << "  QD     : n/a\n";
	}

	fReportView->LoadBuffer(os.str().c_str());
}


inline void HistoGUIUnifiedRoot::PrintRunSummaryToTerminal(const PlotResults& r) const
{
	const double j1v = (std::fabs(r.j1) > 1e-14) ? r.j1 : GetJ1();
	const double j2v = (std::fabs(r.j2) > 1e-14) ? r.j2 : GetJ2();
	const double Eg  = (std::fabs(r.gamma_keV) > 1e-14) ? r.gamma_keV : GetGammaKeV();
	const double sig = (std::fabs(r.sigma) > 1e-14) ? r.sigma : GetSigma();

	std::printf("\n============================================================\n");
	std::printf(" RUN SUMMARY\n");
	std::printf("------------------------------------------------------------\n");
	std::printf(" Inputs:\n");
	std::printf("   E_gamma = %.6f keV\n", Eg);
	std::printf("   j1      = %.6f\n", j1v);
	std::printf("   j2      = %.6f\n", j2v);
	std::printf("   sigma   = %.6f\n", sig);

	std::printf(" Fit results:\n");
	std::printf("   A0 = %.6f\n", A0);
	std::printf("   A2 = %.6f\n", A2E);
	std::printf("   A4 = %.6f\n", A4E);

	if (fHasBestChi2) {
		std::printf(" Chi2 minimum:\n");
		std::printf("   atan(delta)_min = %.6f rad\n", fLastBestAtanDelta);
	} else {
		std::printf(" Chi2 minimum:\n");
		std::printf("   n/a\n");
	}

	if (fHasQD) {
		std::printf(" Detector attenuation:\n");
		std::printf("   QD2 = %.6f\n", fLastQD2);
		std::printf("   QD4 = %.6f\n", fLastQD4);
	} else {
		std::printf(" Detector attenuation:\n");
		std::printf("   QD2/QD4 = n/a\n");
	}
	std::printf("============================================================\n\n");
}

//==============================================================================
// Crosshair helpers
//==============================================================================
inline void HistoGUIUnifiedRoot::ClearCrosshair()
{
	fCrosshairActive = false;
	fCrosshairViewMode = -1;
	fCrossX = 0.0;
	fCrossY = 0.0;

	delete fCrossVLine; fCrossVLine = nullptr;
	delete fCrossHLine; fCrossHLine = nullptr;
}

inline bool HistoGUIUnifiedRoot::GetCurrentPlotUserCoordsFromPixel(Int_t px, Int_t py, double& x, double& y) const
{
	if (!fECanvas || !fECanvas->GetCanvas()) return false;
	TVirtualPad* pad = fECanvas->GetCanvas();
	if (!pad) return false;
	if (!fFrameHist) return false;

	pad->cd();
	x = pad->PadtoX(pad->AbsPixeltoX(px));
	y = pad->PadtoY(pad->AbsPixeltoY(py));
	return true;
}

inline void HistoGUIUnifiedRoot::DrawCrosshairOverlay()
{
	if (!fCrosshairActive) return;
	if (fCrosshairViewMode != fViewMode) return;
	if (!fECanvas || !fECanvas->GetCanvas()) return;
	if (!fFrameHist) return;

	TCanvas* c = fECanvas->GetCanvas();
	c->cd();

	const double xmin = fFrameHist->GetXaxis()->GetXmin();
	const double xmax = fFrameHist->GetXaxis()->GetXmax();
	const double ymin = fFrameHist->GetMinimum();
	const double ymax = fFrameHist->GetMaximum();

	delete fCrossVLine; fCrossVLine = nullptr;
	delete fCrossHLine; fCrossHLine = nullptr;

	fCrossVLine = new TLine(fCrossX, ymin, fCrossX, ymax);
	fCrossHLine = new TLine(xmin, fCrossY, xmax, fCrossY);

	fCrossVLine->SetLineColor(kBlue + 2);
	fCrossHLine->SetLineColor(kBlue + 2);
	fCrossVLine->SetLineStyle(2);
	fCrossHLine->SetLineStyle(2);
	fCrossVLine->SetLineWidth(2);
	fCrossHLine->SetLineWidth(2);

	fCrossVLine->Draw();
	fCrossHLine->Draw();

	c->Modified();
	c->Update();
}

inline void HistoGUIUnifiedRoot::HandleCanvasEvent(int event, int px, int py, TObject* /*selected*/)
{
	if (event != kButton1Down) return;
	if (!fFrameHist) return;

	double x = 0.0, y = 0.0;
	if (!GetCurrentPlotUserCoordsFromPixel(px, py, x, y)) return;

	fCrosshairActive   = true;
	fCrosshairViewMode = fViewMode;
	fCrossX            = x;
	fCrossY            = y;

	if (fViewMode == kViewAngular) {
		std::printf("Canvas click [Angular]: theta = %.6f deg, y = %.6f\n", x, y);
	} else {
		std::printf("Canvas click [Chi2]: x = %.6f, y = %.6f\n", x, y);
	}

	std::ostringstream os;
	os << std::fixed << std::setprecision(6);
	if (fViewMode == kViewAngular) {
		os << "Crosshair: theta=" << x << " deg, y=" << y;
	} else {
		os << "Crosshair: x=" << x << ", y=" << y;
	}
	SetStatus(os.str());

	DrawCurrentView(false);
}


//==============================================================================
// ROOT GUI event handling
//==============================================================================
inline Bool_t HistoGUIUnifiedRoot::ProcessMessage(Longptr_t msg, Longptr_t parm1, Longptr_t /*parm2*/)
{
	switch (GET_MSG(msg)) {
		case kC_COMMAND:
			switch (GET_SUBMSG(msg)) {
				case kCM_BUTTON:
					switch ((int)parm1) {
						case kID_BrowseData: {
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
										     SetDataFileName(fi.fFilename);
									     }
									     if (fi.fIniDir) {
										     dir = fi.fIniDir;
									     }
									     break;
								     }

						case kID_RunCompute: {
									     const std::string fn = GetDataFileName();
									     if (fn.empty()) {
										     SetStatus("Select an angular data file first.");
										     break;
									     }
									     if (gSystem && gSystem->AccessPathName(fn.c_str()) != 0) {
										     SetStatus("File not found / not accessible.");
										     break;
									     }

									     fPendingReq.j1 = GetJ1();
									     fPendingReq.j2 = GetJ2();
									     fPendingReq.gamma_keV = GetGammaKeV();
									     fPendingReq.sigma = GetSigma();
									     fPendingReq.dataFile = fn;

									     fRunRequested = true;
									     ClearCrosshair();
									     SetStatus("Run requested... computing in AD3.cxx");
									     break;
								     }

						case kID_Mode1: {
									SetViewMode(kViewAngular);
									break;
								}

						case kID_Mode2: {
									SetViewMode(kViewChi2);
									break;
								}

						case kID_Mode3Close: {
									     Close();
									     break;
								     }

						case kID_Redraw: {
									 Redraw();
									 break;
								 }

						case kID_Reset: {
									ResetView();
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

#endif // GUI_UNIFIED_ROOT_H
