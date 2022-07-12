
#include "global.h"
//To compile : g++ AD.cxx -o {Input Executable Name} -lX11

#include "Coeff.h"

using namespace std;



double CG(double j1, double j2, double J, double m1, double m2, double M){
	//recall that j1,m1 + j2,m2 = J,M

//	printf("-----------------\n");
	if(M != m1 + m2) return 0;

	double Jmin  =  abs(j1 - j2);
	double Jmax  =  j1+j2;

	
	if(J < Jmin || Jmax < J) return 0;
	


	double a0  = (2*J+1.0)*tgamma(J+j1-j2+1) * tgamma(J-j1+j2+1) * tgamma(j1+j2-J+1)/tgamma(J+j1+j2+1.0 +1);
	double A0 = sqrt(a0);


	double a = tgamma(J+M+1)   *tgamma(J-M+1);
	double a1= tgamma(j1+m1+1) *tgamma(j1-m1+1);
	double a2= tgamma(j2+m2+1) *tgamma(j2-m2+1);

	double A = sqrt( a  * a1 * a2);


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
	return A0*A*cg;
}


double ThreeJsym(double j1, double m1, double j2, double m2, double j3, double m3){
	//[j1 j2 j3] = (-1)^(j1-j2-m3)/ sqrt(2*j3+1) * CG(j3, -m3, j1, m1, j2, m2)
	//[m1 m2 m3]
	
	
											  //J,M,j1,m1,j2,m2)	
	//return pow(-1,j1 -j2 -m3)/sqrt(2*j3+1)*CG(j3,-m3,j1,m1,j2,m2);
	double threej =  pow(-1,j1 -j2 -m3)/sqrt(2*j3+1)*CG(j3,-m3,j1,m1,j2,m2);
	
	//double threej =  pow(-1,j1 -j2 -m3)/sqrt(2*j3+1)*CG(j1,j2,j3,m1,m2,-m3);
	printf("---------\n");
	printf("3J= %lf\n",threej);
	
	printf("---------\n");

	return threej;
}

double SixJsym(double j1, double j2, double j3, double j4, double j5, double j6){

//--------------------------------------------------------------------------------//
	// The six j symbol describes the coupling between j1 j2 and j3 to J - j1.
	// essentially a triangle of angular momentum rules between these. 
	// j1 = j1
	// j2 = j2
	// j3 = j1 + j2
	// j4 = j3
	// j5 = J = j1 + j2 + j3
	// j6 = j2 + j3
// ----------------------------------------------------------------------------- //
// the following conditions check the triangle selection rules

	if( j3 < abs(j1 - j2) || j1 + j2 < j3) return 0; 
	if( j6 < abs(j2 - j4) || j2 + j4 < j6) return 0; 
	if( j5 < abs(j1 - j6) || j1 + j6 < j5) return 0; 
	if( j5 < abs(j3 - j4) || j3 + j4 < j5) return 0; 
// now that they have been checked, we can go ahead and calculate sixJ.

	double sixj = 0.0;
	float m1 = -j1;
	float m2 = -j2;
	float m3 = -j3;
	float m4 = -j4;
	float m5 = -j5;
	float m6 = -j6;
printf("here\n");
	for(; m1 <= j1; m1 = m1 +1){
		for(; m2 <= j2; m2 = m2 +1){
			for(; m3 <= j3; m3 = m3 + 1){
				for(; m4 <= j4; m4 = m4 +1){
					for(; m5 <= j5; m5 = m5 + 1){
						for(; m6 <= j6; m6 = m6 +1){
	
							double h = (j1 - m1) + (j2 - m2) + (j3 -m3) + (j4 - m4) + (j5 - m5) + (j6 - m6);
							double b1 = ThreeJsym(j1, -m1, j2, -m2, j3, -m3);
							double b2 = ThreeJsym(j1, m1, j5, -m5, j6, m6); 
							double b3 = ThreeJsym(j4, m4, j2, m2,  j6, -m6);
							double b4 = ThreeJsym(j4, -m4, j5, m5, j3, m3);					
							double b = b1 * b2 * b3 * b4;
		

							printf("h = %lf\n", h);						
							printf("b1 = %lf\n", b1);						
							printf("b2 = %lf\n", b2);						
							printf("b3 = %lf\n", b3);						
							printf("b4 = %lf\n", b4);						
							printf("b = %lf\n", b);						

							sixj += pow(-1,h)*b;
	
						}
					}
				}
			}
		}
	}

	return sixj;
}


int main(int argc, char ** argv){

	//if mod (2*J1,2) = 1 do this
//	double A0[6] = {1,1,2,0,0,0};
//	double A1[6] = {1,1,4,0,0,0};
	
	//if mod(2*j1,20 = 0 do this
	
	double j1 = 1.;
	double j2 = 2.;


	double A0[6] = {j1,j1,2,.5,-.5,0};
	double A1[6] = {j1,j1,2,-.5,.5,0};
	double A2[6] = {j1,j1,4,-.5,.5,0};
	double A3[6] = {j1,j1,4,.5,-.5,0};

	double cgr0= CG(j1,j1,2,.5,-.5,0);
	double cgr1= CG(j1,j1,2,-.5,.5,0);
	double cgr2= CG(j1,j1,4,.5,-.5,0);
	double cgr3= CG(j1,j1,4,-.5,.5,0);
	

	printf("------\n");
//	printf("CGR0 = %lf\n",cgr0);
//	printf("CGR1 = %lf\n",cgr1);
//	printf("CGR2 = %lf\n",cgr2);
//	printf("CGR3 = %lf\n",cgr3);

	double sixj = SixJsym(j1,j1,1.,1.,2.,2.);
	double sixj2 = SixJsym(j1,j1,1.,2.,2.,2.);
	double sixj3 = SixJsym(j1,j1,2.,2.,2.,2.); 



	printf("6-J = %lf\n", sixj);



	return 0; 
}



