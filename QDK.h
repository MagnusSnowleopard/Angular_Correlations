//#####################################################################
//
//This file is a geometric calculation of the efficiency of radiation
//in germanium crystals as a function of energy, and the crystal
//geometry of the experiment.
//
//#####################################################################

#include "global.h"

using namespace std;

namespace {

static inline void ComputeQDKValues(double Energy,
		double radius,
		double distance,
		double thickness,
		double& Tau,
		double& QD2,
		double& QD4)
{
	const double E_mev = Energy / 1000.0;
	const double E_log = log(E_mev);

	const double EL1 = E_log;
	const double EL2 = EL1 * EL1;
	const double EL3 = EL1 * EL2;
	const double EL4 = EL2 * EL2;
	const double EL5 = EL4 * EL1;

	const double TT = -1.1907 - 0.5372*EL1 - 0.0438*EL2 + 0.0218*EL3 + 0.0765*EL4 + 0.0095*EL5;
	Tau = exp(TT);

	const double Z1 = radius / (distance + thickness);
	const double Z2 = radius / distance;

	const double alpha = atan(Z1);
	const double gamma = atan(Z2);

	const int loop_length = 1000;

	// Region 1: beta in [0, alpha]
	const double BL = 0.0;
	const double BU = alpha;
	const double delx1 = (BU - BL) / loop_length;

	double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

	for (int i = 0; i <= loop_length; ++i) {
		double A = 1.0;
		if (i != 0 && i != loop_length) A = (i % 2 == 0) ? 2.0 : 4.0;
		const double beta = BL + i * delx1;

		const double cosb = cos(beta);
		const double sinb = sin(beta);
		const double secb = 1.0 / cosb;
		const double c2 = cosb * cosb;
		const double c4 = c2 * c2;
		const double ex1 = exp(-Tau * thickness * secb);

		sum1 += 0.5 * (3.0*c2 - 1.0) * (1.0 - ex1) * sinb * A * delx1;
		sum2 += 0.125 * (35.0*c4 - 30.0*c2 + 3.0) * (1.0 - ex1) * sinb * A * delx1;
		sum3 += (1.0 - ex1) * sinb * A * delx1;
	}

	const double ans1 = sum1 / 3.0;
	const double ans2 = sum2 / 3.0;
	const double ans3 = sum3 / 3.0;

	// Region 2: beta in [alpha, gamma]
	const double LB = alpha;
	const double UB = gamma;
	const double delx2 = (UB - LB) / loop_length;

	double sum4 = 0.0, sum5 = 0.0, sum6 = 0.0;

	for (int i = 0; i <= loop_length; ++i) {
		double A = 1.0;
		if (i != 0 && i != loop_length) A = (i % 2 == 0) ? 2.0 : 4.0;
		const double beta = LB + i * delx2;

		const double cosb = cos(beta);
		const double sinb = sin(beta);
		const double secb = 1.0 / cosb;
		const double cscb = 1.0 / sinb;
		const double c2 = cosb * cosb;
		const double c4 = c2 * c2;
		const double ex2 = exp(-Tau * (radius * cscb - distance * secb));

		sum4 += 0.5 * (3.0*c2 - 1.0) * (1.0 - ex2) * sinb * A * delx2;
		sum5 += 0.125 * (35.0*c4 - 30.0*c2 + 3.0) * (1.0 - ex2) * sinb * A * delx2;
		sum6 += (1.0 - ex2) * sinb * A * delx2;
	}

	const double ans4 = sum4 / 3.0;
	const double ans5 = sum5 / 3.0;
	const double ans6 = sum6 / 3.0;

	const double denom = ans3 + ans6;
	QD2 = (ans1 + ans4) / denom;
	QD4 = (ans2 + ans5) / denom;
}

} // namespace


double QK2(double Energy, double radius, double distance, double thickness){
	double Tau = 0.0, QD2 = 0.0, QD4 = 0.0;
	ComputeQDKValues(Energy, radius, distance, thickness, Tau, QD2, QD4);

	//Now output a file that contains R, D , T , gamma energy, attentuation coeff, q2 and q4
	ofstream fileo;
	fileo.open ("ad.txt");
	fileo << "Radius = " << radius <<" [cm]\n";
	fileo << "Distance = " << distance <<" [cm]\n";
	fileo << "Thickness = " << thickness <<" [cm]\n";
	fileo << "Atten.C = " << Tau <<" [cm^-1]\n";
	fileo << "Gamma_E = " << Energy <<" [KeV]\n";
	fileo << "QD2 = " << QD2 << "\n";
	fileo << "QD4 = " << QD4 << "\n";
	fileo.close();

	return QD2;
}


double QK4(double Energy, double radius, double distance, double thickness){
	double Tau = 0.0, QD2 = 0.0, QD4 = 0.0;
	ComputeQDKValues(Energy, radius, distance, thickness, Tau, QD2, QD4);
	return QD4;
}
