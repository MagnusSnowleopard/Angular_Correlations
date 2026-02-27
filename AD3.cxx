//#################################################################
//  Angular Distribution + Chi^2 minimization
//  Racah coefficients computed via Racah3.h 
//
//  Compile:
//    g++ -std=c++17 AD3.cxx `root-config --cflags --glibs` -lGui -Wno-ignored-attributes -o AD
//
//#################################################################

#include "global.h"
#include "GUI_UNIFIED_ROOT.h"
#include "Racah3.h"

#include <TApplication.h>
#include <TSystem.h>

#include "Functlib.h"
#include "QDK.h"
#include "Cleb.h"


#include <iomanip>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// --------------------------
// Flexible tokenize: comma + whitespace
// --------------------------
static std::vector<std::string> TokenizeFlexible(const std::string& line)
{
	std::string s = line;
	for (size_t i = 0; i < s.size(); ++i) {
		if (s[i] == ',') s[i] = ' ';
	}
	std::istringstream iss(s);
	std::vector<std::string> out;
	std::string tok;
	while (iss >> tok) out.push_back(tok);
	return out;
}

static bool Parse3D(const std::vector<std::string>& row, double& a, double& b, double& c)
{
	if (row.size() < 3) return false;
	char* e1=nullptr; char* e2=nullptr; char* e3=nullptr;
	a = std::strtod(row[0].c_str(), &e1);
	b = std::strtod(row[1].c_str(), &e2);
	c = std::strtod(row[2].c_str(), &e3);
	if (e1 == row[0].c_str() || e2 == row[1].c_str() || e3 == row[2].c_str()) return false;
	return true;
}

// --------------------------
// Load experimental angular data file:
// expected columns: angle(deg)  y  ey
// Returns angles in radians.
// --------------------------
static bool LoadAngularExperimentalFile(const std::string& filename,
		std::vector<double>& theta_rad,
		std::vector<double>& y,
		std::vector<double>& ey,
		std::string& err)
{
	theta_rad.clear();
	y.clear();
	ey.clear();

	std::ifstream f(filename.c_str());
	if (!f.is_open()) {
		err = "Could not open data file: " + filename;
		return false;
	}

	std::string line;
	while (std::getline(f, line)) {
		if (line.empty()) continue;

		// allow comments + skip blank
		std::string t = line;
		t.erase(t.begin(), std::find_if(t.begin(), t.end(),
					[](unsigned char ch){ return !std::isspace(ch); }));
		if (t.empty()) continue;
		if (t[0] == '#') continue;

		std::vector<std::string> row = TokenizeFlexible(line);

		double a=0.0, b=0.0, c=0.0;
		if (!Parse3D(row, a, b, c)) continue; // skips headers/non-numeric rows

		// File angles are degrees (matches your legacy AD.cxx)
		double th = a * M_PI / 180.0;

		theta_rad.push_back(th);
		y.push_back(b);
		ey.push_back(c);
	}

	if (theta_rad.empty() || theta_rad.size() != y.size() || y.size() != ey.size()) {
		err = "No valid numeric rows found (need 3 columns: angle(deg) y ey).";
		return false;
	}

	return true;
}

// --------------------------
// Solve 3x3 normal equations for Legendre fit coefficients A0,A2,A4
// (unweighted, matches your original)
// --------------------------
static bool FitLegendreA0A2A4(const std::vector<double>& theta_rad,
		const std::vector<double>& y,
		double& A0, double& A2, double& A4,
		std::string& err)
{
	if (theta_rad.size() < 3 || theta_rad.size() != y.size()) {
		err = "Not enough angular points for A0/A2/A4 fit.";
		return false;
	}

	double a1=0.,a2=0.,a3=0.,a4s=0.,a5=0.,a6=0.,a7=0.,a8=0.,a9=0.,a10=0.,a11=0.,a12=0.;

	for (size_t i = 0; i < theta_rad.size(); ++i) {
		const double th = theta_rad[i];
		const double c  = std::cos(th);
		const double c2 = c*c;
		const double c4 = c2*c2;

		const double P2 = (1.5*c2 - 0.5);
		const double P4 = (35.0/8.0*c4 - 30.0/8.0*c2 + 3.0/8.0);

		a1  += 1.0;
		a2  += P2;
		a3  += P4;
		a4s += y[i];

		a5  += P2;
		a6  += P2*P2;
		a7  += P2*P4;
		a8  += P2*y[i];

		a9  += P4;
		a10 += P2*P4;
		a11 += P4*P4;
		a12 += P4*y[i];
	}

	double mat[3][4] = {
		{a1,  a2,  a3,  a4s},
		{a5,  a6,  a7,  a8 },
		{a9,  a10, a11, a12}
	};

	// partial pivot
	for (int i = 0; i < 3; ++i) {
		int piv = i;
		for (int r = i+1; r < 3; ++r) {
			if (std::fabs(mat[r][i]) > std::fabs(mat[piv][i])) piv = r;
		}
		if (piv != i) {
			for (int k = 0; k < 4; ++k) std::swap(mat[i][k], mat[piv][k]);
		}
		if (std::fabs(mat[i][i]) < 1e-18) {
			err = "Singular normal-equation matrix in Legendre fit (check input angles).";
			return false;
		}
	}

	// elimination
	for (int i = 0; i < 2; ++i) {
		for (int r = i+1; r < 3; ++r) {
			double f = mat[r][i] / mat[i][i];
			for (int k = i; k < 4; ++k) mat[r][k] -= f*mat[i][k];
		}
	}

	// back-sub
	double x[3] = {0.,0.,0.};
	for (int i = 2; i >= 0; --i) {
		double s = mat[i][3];
		for (int j = i+1; j < 3; ++j) s -= mat[i][j]*x[j];
		x[i] = s / mat[i][i];
	}

	A0 = x[0];
	A2 = x[1];
	A4 = x[2];

	if (std::fabs(A0) < 1e-18) {
		err = "Fit returned A0 ~ 0 (cannot normalize).";
		return false;
	}

	return true;
}

// --------------------------
// Compute Bk11, Bk12 from sigma and j1 (legacy logic)
// IMPORTANT: uses integer parity, not pow(-1, double)
// --------------------------
static void ComputeBk(double j1, double sigma, double& Bk11, double& Bk12)
{
	Bk11 = 0.0;
	Bk12 = 0.0;

	const double js = std::sqrt(2*j1 + 1);
	const double j12 = 2*j1;
	const int IS = ((int)std::llround(j12)) % 2;

	if (sigma == 0.0) {
		if (IS == 1) {
			double A[6] = {j1, j1, 2.0, 0.5, -0.5, 0.0};
			double cg12 = CG2(A);

			A[3] = -0.5; A[4] = 0.5;
			double cg22 = CG2(A);

			A[2] = 4.0;
			double cg24 = CG2(A);

			A[3] = 0.5; A[4] = -0.5;
			double cg14 = CG2(A);

			// legacy placeholders
			long long I1i = 0;
			long long I2i = 0;
			double ph1 = (I1i % 2 == 0) ? 1.0 : -1.0;
			double ph2 = (I2i % 2 == 0) ? 1.0 : -1.0;

			Bk11 = 0.5*js*(ph1*cg12 + ph2*cg22);
			Bk12 = 0.5*js*(ph1*cg14 + ph2*cg24);

		} else {
			double A[6] = {j1, j1, 2.0, 0.0, 0.0, 0.0};
			double cg1 = CG2(A);
			// pow(-1,j1) is unsafe for half-integers; approximate parity via 2j
			long long tw = (long long)std::llround(2.0*j1);
			// (-1)^(j1) not well-defined for half integer; legacy expression used it anyway.
			// Keep phase = +1 for even 2j, -1 for odd 2j (i.e., (-1)^(2j1/2) -> (-1)^j1)
			double ph = ((tw/2) % 2 == 0) ? 1.0 : -1.0;
			Bk11 = ph * js * cg1;

			A[2] = 4.0;
			double cg2 = CG2(A);
			Bk12 = ph * js * cg2;
		}
		return;
	}

	// sigma != 0 : gaussian W(m)
	const double sigsq = sigma*sigma;
	const double j14 = 4*j1;

	double sum1 = 0.0;
	for (int i = 0; i <= (int)std::llround(j14); i += 2) {
		double m1 = 0.5*(i - j12);
		double ex = std::exp(-(m1*m1)/(2*sigsq));
		sum1 += ex;
	}
	const double cn1 = (sum1 != 0.0) ? (1.0/sum1) : 0.0;

	const double sfact = std::sqrt(2*j1 + 1);

	for (int K = 2; K <= 4; K += 2) {
		for (int m = 0; m <= (int)std::llround(j14); m += 2) {
			double m1 = 0.5*(m - j12);
			double ex = std::exp(-(m1*m1)/(2*sigsq));

			double II = j1 - m1;
			long long IIi = std::llround(II);
			double phase = (IIi % 2 == 0) ? 1.0 : -1.0;

			double cgg = CGcoeff(K, 0, j1, m1, j1, -m1);

			double tTerm = cn1 * ex * phase * sfact * cgg;

			if (K == 2) Bk11 += tTerm;
			if (K == 4) Bk12 += tTerm;
		}
	}
}

// --------------------------
// Main computation engine called from GUI Run
// Fills PlotResults for GUI display
// Uses Racah.h ComputeRK(J1,J2) instead of RkTable.csv.
// --------------------------
static bool RunADComputation(const HistoGUIUnifiedRoot::RunRequest& req,
		double detradius,
		double targetdistance,
		double detthickness,
		HistoGUIUnifiedRoot::PlotResults& out,
		std::string& err)
{
	// Validate
	if (req.dataFile.empty()) {
		err = "No data file selected.";
		return false;
	}

	// Load experimental data
	std::vector<double> theta_rad, yraw, eyraw;
	if (!LoadAngularExperimentalFile(req.dataFile, theta_rad, yraw, eyraw, err)) {
		return false;
	}

	// Fit Legendre to experimental data -> A0,A2,A4
	double A0E=1.0, A2E=0.0, A4E=0.0;
	if (!FitLegendreA0A2A4(theta_rad, yraw, A0E, A2E, A4E, err)) {
		return false;
	}

	// Detector attenuation coefficients using GUI gamma energy
	const double Energy = req.gamma_keV;
	if (!(Energy > 0.0)) {
		err = "Gamma energy must be > 0.";
		return false;
	}

	const double QD2 = QK2(Energy, detradius, targetdistance, detthickness);
	const double QD4 = QK4(Energy, detradius, targetdistance, detthickness);
	/*	std::cout << std::fixed << std::setprecision(6)
		<< "[QDK] E_keV=" << Energy
		<< "  R=" << detradius
		<< "  D=" << targetdistance
		<< "  T=" << detthickness
		<< "  => QD2=" << QD2
		<< "  QD4=" << QD4 << "\n";
		*/	// --- Racah coefficients computed from Racah.h ---
	double rk01=0, rk11=0, rk21=0, rk02=0, rk12=0, rk22=0;
	int Lbase = -1;
	try {
		RACAH_COEFFS::RKSet r = RACAH_COEFFS::ComputeRK(req.j1, req.j2);
		rk01 = r.rk01; rk11 = r.rk11; rk21 = r.rk21;
		rk02 = r.rk02; rk12 = r.rk12; rk22 = r.rk22;
		Lbase = r.Lbase;
		(void)Lbase;
	} catch (const std::exception& e) {
		err = std::string("Racah ComputeRK failed: ") + e.what();
		return false;
	}

	// Compute Bk11,Bk12 from sigma
	double Bk11=0.0, Bk12=0.0;
	if (req.sigma < 0.0) {
		err = "Sigma must be >= 0.";
		return false;
	}
	ComputeBk(req.j1, req.sigma, Bk11, Bk12);

	// Experimental normalized coefficients
	const double a2E = A2E / A0E;
	const double a4E = A4E / A0E;
	/*	std::cout << std::fixed << std::setprecision(6)
		<< "[Legendre Fit] A0=" << A0E
		<< "  A2=" << A2E
		<< "  A4=" << A4E
		<< "   (a2=A2/A0=" << (A2E / A0E)
		<< ", a4=A4/A0=" << (A4E / A0E) << ")\n";
		*/	// Chi2 scan setup
	const double delta_min = -M_PI/2.0;
	const double delta_max =  M_PI/2.0;
	const double step = 0.01;

	const int points = (int)((delta_max - delta_min)/step);

	std::vector<double> chisqr;
	std::vector<double> tdelta;
	chisqr.reserve(points);
	tdelta.reserve(points);

	// Chi2 scan comparing theory vs experimental-fit curve (legacy style)
	for (int i = 0; i < points; ++i) {
		double atan_delta = i*step + delta_min;
		double delta = std::tan(atan_delta);

		double X2_total = 0.0;

		for (size_t j = 0; j < theta_rad.size(); ++j) {
			const double Tangle = theta_rad[j];
			const double Y_err  = eyraw[j];

			const double denom = (theta_rad.size() > 1 ? (theta_rad.size()-1) : 1) * (Y_err*Y_err);
			if (denom <= 0.0) continue;

			const double rd2T = (rk01 + 2*delta*rk11 + delta*delta*rk21) / (1 + delta*delta);
			const double rd4T = (rk02 + 2*delta*rk12 + delta*delta*rk22) / (1 + delta*delta);

			const double c  = std::cos(Tangle);
			const double c2 = c*c;
			const double c4 = c2*c2;

			const double P2 = (1.5*c2 - 0.5);
			const double P4 = (35.0/8.0*c4 - 30.0/8.0*c2 + 3.0/8.0);

			const double YT0 = 1.0;
			const double YT2 = QD2 * Bk11 * rd2T * P2;
			const double YT4 = QD4 * Bk12 * rd4T * P4;
			const double YT  = YT0 + YT2 + YT4;

			const double YE  = 1.0 + a2E * P2 + a4E * P4;

			X2_total += ((YT - YE)*(YT - YE)) / denom;
		}

		if (X2_total <= 0.0) X2_total = 1e-30; // avoid log(0)
		chisqr.push_back(std::log(X2_total));
		tdelta.push_back(atan_delta);
	}

	// Prepare normalized angular y for display
	std::vector<double> ynorm;
	ynorm.reserve(yraw.size());
	for (size_t i = 0; i < yraw.size(); ++i) ynorm.push_back(yraw[i] / A0E);

	// Fill outputs for GUI
	out = HistoGUIUnifiedRoot::PlotResults{};
	out.dangler = theta_rad;
	out.dydatas = ynorm;
	out.deydata = eyraw;

	out.fitA0 = A0E;
	out.fitA2 = A2E;
	out.fitA4 = A4E;
	out.hasFit = true;

	out.tdelta = tdelta;
	out.chisqr = chisqr;

	out.j1 = req.j1;
	out.j2 = req.j2;
	out.gamma_keV = req.gamma_keV;
	out.sigma = req.sigma;

	return true;
}

#ifndef __CINT__

int main(int argc, char** argv)
{
	TApplication app("AD_gui_app", &argc, argv);

	// ------------------------------------------------------------
	// Detector geometry stays in AD.cxx (your request)
	// ------------------------------------------------------------
	double detradius      = 3.0;
	double targetdistance = 4.0;
	double detthickness   = 5.0;

	// GUI defaults
	double default_j1       = 7.5;
	double default_j2       = 5.5;
	double default_gammaKeV = 426.0;
	double default_sigma    = 1.0;

	// ------------------------------------------------------------
	// Create GUI
	// ------------------------------------------------------------
	HistoGUIUnifiedRoot* gui = new HistoGUIUnifiedRoot(gClient->GetRoot(), 1400, 920);

	gui->SetDetectorGeometry(detradius, targetdistance, detthickness);

	gui->SetJ1(default_j1);
	gui->SetJ2(default_j2);
	gui->SetGammaKeV(default_gammaKeV);
	gui->SetSigma(default_sigma);

	gui->Init();
	gui->SetViewMode(HistoGUIUnifiedRoot::kViewAngular);

	// ------------------------------------------------------------
	// GUI-driven main loop (no terminal menu)
	// ------------------------------------------------------------
	while (!gui->ShouldQuitProgram()) {
		gSystem->ProcessEvents();
		gui->PollCanvasClick();
		gSystem->Sleep(20);

		HistoGUIUnifiedRoot::RunRequest req;
		if (!gui->ConsumeRunRequest(req)) continue;

		gui->SetStatus("Computing...");

		HistoGUIUnifiedRoot::PlotResults results;
		std::string err;

		bool ok = RunADComputation(req, detradius, targetdistance, detthickness, results, err);

		if (!ok) {
			gui->SetStatus(err.empty() ? "Computation failed." : err);
			continue;
		}

		gui->LoadResults(results);
		gui->SetStatus("Done. Use 1/2 to switch views; 3 to close.");
	}

	delete gui;
	gui = nullptr;

	return 0;
}

#endif
