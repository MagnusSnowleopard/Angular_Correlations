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
#include "QDK2.h"
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
		const std::vector<double>& ey,
		double& A0, double& A2, double& A4,
		double& a2, double& a4,
		double& a2Err, double& a4Err,
		std::string& err)
{
	if (theta_rad.size() < 3 || theta_rad.size() != y.size() || y.size() != ey.size()) {
		err = "Not enough angular points for A0/A2/A4 fit.";
		return false;
	}

	double M[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
	double b[3] = {0.,0.,0.};

	for (size_t i = 0; i < theta_rad.size(); ++i) {
		const double th = theta_rad[i];
		const double c  = std::cos(th);
		const double c2 = c*c;
		const double c4 = c2*c2;

		const double P2 = (1.5*c2 - 0.5);
		const double P4 = (35.0/8.0*c4 - 30.0/8.0*c2 + 3.0/8.0);
		const double f[3] = {1.0, P2, P4};

		double w = 1.0;
		if (ey[i] > 0.0 && std::isfinite(ey[i])) w = 1.0/(ey[i]*ey[i]);

		for (int r = 0; r < 3; ++r) {
			b[r] += w * f[r] * y[i];
			for (int cidx = 0; cidx < 3; ++cidx) M[r][cidx] += w * f[r] * f[cidx];
		}
	}

	double Aaug[3][6] = {
		{M[0][0], M[0][1], M[0][2], 1.0, 0.0, 0.0},
		{M[1][0], M[1][1], M[1][2], 0.0, 1.0, 0.0},
		{M[2][0], M[2][1], M[2][2], 0.0, 0.0, 1.0}
	};

	for (int i = 0; i < 3; ++i) {
		int piv = i;
		for (int r = i+1; r < 3; ++r) if (std::fabs(Aaug[r][i]) > std::fabs(Aaug[piv][i])) piv = r;
		if (piv != i) for (int k = 0; k < 6; ++k) std::swap(Aaug[i][k], Aaug[piv][k]);
		if (std::fabs(Aaug[i][i]) < 1e-18) {
			err = "Singular normal-equation matrix in Legendre fit (check input angles).";
			return false;
		}
		double diag = Aaug[i][i];
		for (int k = 0; k < 6; ++k) Aaug[i][k] /= diag;
		for (int r = 0; r < 3; ++r) {
			if (r == i) continue;
			double f = Aaug[r][i];
			for (int k = 0; k < 6; ++k) Aaug[r][k] -= f * Aaug[i][k];
		}
	}

	double Minv[3][3] = {
		{Aaug[0][3], Aaug[0][4], Aaug[0][5]},
		{Aaug[1][3], Aaug[1][4], Aaug[1][5]},
		{Aaug[2][3], Aaug[2][4], Aaug[2][5]}
	};

	A0 = Minv[0][0]*b[0] + Minv[0][1]*b[1] + Minv[0][2]*b[2];
	A2 = Minv[1][0]*b[0] + Minv[1][1]*b[1] + Minv[1][2]*b[2];
	A4 = Minv[2][0]*b[0] + Minv[2][1]*b[1] + Minv[2][2]*b[2];

	if (std::fabs(A0) < 1e-18) {
		err = "Fit returned A0 ~ 0 (cannot normalize).";
		return false;
	}

	a2 = A2 / A0;
	a4 = A4 / A0;

	const double varA0 = Minv[0][0];
	const double varA2 = Minv[1][1];
	const double varA4 = Minv[2][2];
	const double covA0A2 = Minv[0][1];
	const double covA0A4 = Minv[0][2];

	const double A0sq = A0*A0;
	const double A0cu = A0sq*A0;
	const double var_a2 = varA2/A0sq + (A2*A2)*varA0/(A0sq*A0sq) - 2.0*A2*covA0A2/A0cu;
	const double var_a4 = varA4/A0sq + (A4*A4)*varA0/(A0sq*A0sq) - 2.0*A4*covA0A4/A0cu;

	a2Err = (var_a2 > 0.0 && std::isfinite(var_a2)) ? std::sqrt(var_a2) : 0.0;
	a4Err = (var_a4 > 0.0 && std::isfinite(var_a4)) ? std::sqrt(var_a4) : 0.0;

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

	// Fit Legendre to experimental data -> A0,A2,A4 and normalized a2,a4
	double A0E=1.0, A2E=0.0, A4E=0.0;
	double a2E=0.0, a4E=0.0, a2Err=0.0, a4Err=0.0;
	if (!FitLegendreA0A2A4(theta_rad, yraw, eyraw, A0E, A2E, A4E, a2E, a4E, a2Err, a4Err, err)) {
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

	// Experimental normalized coefficients were computed in the fit stage
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
	std::vector<double> chi2raw;
	std::vector<double> tdelta;
	chisqr.reserve(points);
	chi2raw.reserve(points);
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
		chi2raw.push_back(X2_total);
		chisqr.push_back(std::log(X2_total));
		tdelta.push_back(atan_delta);
	}

	// Estimate delta uncertainty from chi2_min + 1 crossing
	double bestDelta = 0.0;
	double bestDeltaErr = 0.0;
	bool hasDeltaErr = false;
	if (!chi2raw.empty() && chi2raw.size() == tdelta.size()) {
		size_t imin = 0;
		for (size_t i = 1; i < chi2raw.size(); ++i) if (chi2raw[i] < chi2raw[imin]) imin = i;
		bestDelta = std::tan(tdelta[imin]);
		const double target = chi2raw[imin] + 1.0;

		int ileft = (int)imin;
		while (ileft > 0 && chi2raw[(size_t)ileft] <= target) --ileft;
		int iright = (int)imin;
		while (iright + 1 < (int)chi2raw.size() && chi2raw[(size_t)iright] <= target) ++iright;

		if (ileft < (int)imin && iright > (int)imin && ileft >= 0 && iright < (int)tdelta.size()) {
			const double dleft = std::tan(tdelta[(size_t)ileft]);
			const double dright = std::tan(tdelta[(size_t)iright]);
			bestDeltaErr = 0.5 * std::fabs(dright - dleft);
			hasDeltaErr = std::isfinite(bestDeltaErr) && bestDeltaErr > 0.0;
		}
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
	out.a2 = a2E;
	out.a4 = a4E;
	out.a2Err = a2Err;
	out.a4Err = a4Err;
	out.hasFit = true;
	out.hasCoeffErrors = std::isfinite(a2Err) && std::isfinite(a4Err) && (a2Err > 0.0 || a4Err > 0.0);

	out.tdelta = tdelta;
	out.chisqr = chisqr;

	out.j1 = req.j1;
	out.j2 = req.j2;
	out.gamma_keV = req.gamma_keV;
	out.sigma = req.sigma;
	out.bestDelta = bestDelta;
	out.bestDeltaErr = bestDeltaErr;
	out.hasDeltaError = hasDeltaErr;

	out.qd2 = QD2;
	out.qd4 = QD4; 
	out.hasQD = std::isfinite(QD2) && std::isfinite(QD4);
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
