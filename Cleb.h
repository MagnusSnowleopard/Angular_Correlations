
#include "global.h"
//To compile : g++ AD.cxx -o {Input Executable Name} -lX11

//#include "Coeff.h"

using namespace std;



//double CG(double J, double M, double j1, double m1, double j2, double m2){
double CG(double j1, double j2, double J, double m1, double m2, double M){
	//recall that j1,m1 + j2,m2 = J,M

	if(M != m1 + m2) return 0;

	double Jmin  =  abs(j1 - j2);
	double Jmax  =  j1+j2;

	
	if(J < Jmin || Jmax < J) return 0;
	
	double a0  = (2*J+1.0)*tgamma(J+j1-j2+1) * tgamma(J-j1+j2+1) * tgamma(j1+j2-J+1)/tgamma(J+j1+j2+1.0 +1);
	double A0 = sqrt(a0);


	double a = tgamma(J+M+1)   *tgamma(J-M+1);
	double a1= tgamma(j1+m1+1) *tgamma(j1-m1+1);
	double a2= tgamma(j2+m2+1) *tgamma(j2-m2+1);

	double A1 = sqrt( a  * a1 * a2);
	int pmax = min( min(j1+j2-J,j1-m1),j2 + m2);
	

	double cg = 0.;

	for( int p =0; p<=pmax;p++){


		double p1 = tgamma(j1+j2-J-p+1);
		double p2 = tgamma(j1-m1-p+1);
		double p3 = tgamma(j2+m2-p+1);
		double p4 = tgamma(J -j2 +m1 +p+1);
		double p5 = tgamma(J -j1 -m2 +p+1);
		double t = pow(-1,p)/(tgamma(p+1) * p1 * p2 * p3 * p4 * p5);

		cg += t;
	}
	return A0*A1*cg;
}


double CG2(double A[6]){
	//recall that j1,m1 + j2,m2 = J,M
	double j1 = A[0];
	double j2 = A[1]; 
	double J = A[2];
	double m1 = A[3];
	double m2 = A[4];
	double M = A[5]; 
	if(M != m1 + m2) return 0;

	double Jmin  =  abs(j1 - j2);
	double Jmax  =  j1+j2;

	
	if(J < Jmin || Jmax < J) return 0;
	


	double a0  = (2*J+1.0)*tgamma(J+j1-j2+1) * tgamma(J-j1+j2+1) * tgamma(j1+j2-J+1)/tgamma(J+j1+j2+1.0 +1);
	double A0 = sqrt(a0);


	double a = tgamma(J+M+1)   *tgamma(J-M+1);
	double a1= tgamma(j1+m1+1) *tgamma(j1-m1+1);
	double a2= tgamma(j2+m2+1) *tgamma(j2-m2+1);

	double A1 = sqrt( a  * a1 * a2);
	int pmax = min( min(j1+j2-J,j1-m1),j2 + m2);
	

	double cg = 0.;

	for( int p =0; p<=pmax;p++){

		double p1 = tgamma(j1+j2-J-p+1);
		double p2 = tgamma(j1-m1-p+1);
		double p3 = tgamma(j2+m2-p+1);
		double p4 = tgamma(J -j2 +m1 +p+1);
		double p5 = tgamma(J -j1 -m2 +p+1);
		double t = pow(-1,p)/(tgamma(p+1) * p1 * p2 * p3 * p4 * p5);

		cg += t;
	}
	return A0*A1*cg;
}
/*

int main(int argc, char ** argv){

	//if mod (2*J1,2) = 1 do this
//	double A0[6] = {1,1,2,0,0,0};
//	double A1[6] = {1,1,4,0,0,0};
	
	//if mod(2*j1,20 = 0 do this
	
	double j1 = 2.;
	double j2 = 2.;


	double A0[6] = {j1,j1,2,.5,-.5,0};
	double A1[6] = {j1,j1,2,-.5,.5,0};
	double A2[6] = {j1,j1,4,-.5,.5,0};
	double A3[6] = {j1,j1,4,.5,-.5,0};

	double cgr0= CG(j1,j1,2,.5,-.5,0);
	double cgr1= CG(j1,j1,2,-.5,.5,0);
	double cgr2= CG(j1,j1,4,.5,-.5,0);
	double cgr3= CG(j1,j1,4,-.5,.5,0);
	

	printf("CGR0 = %lf\n",cgr0);
	printf("CGR1 = %lf\n",cgr1);
	printf("CGR2 = %lf\n",cgr2);
	printf("CGR3 = %lf\n",cgr3);




	return 0; 
}

*/

