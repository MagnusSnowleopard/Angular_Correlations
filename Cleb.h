
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



double ThreeJSymbol(double J1, double m1, double J2, double m2, double J3, double m3){

  // ( J1 J2 J3 ) = (-1)^(J1-J2 - m3)/ sqrt(2*J3+1) * CGcoeff(J3, -m3, J1, m1, J2, m2)
  // ( m1 m2 m3 )

  return pow(-1, J1 - J2 - m3)/sqrt(2*J3+1) * CG(J3, -m3, J1, m1, J2, m2);

}

double TJ2(double A[6]){

  // ( J1 J2 J3 ) = (-1)^(J1-J2 - m3)/ sqrt(2*J3+1) * CGcoeff(J3, -m3, J1, m1, J2, m2)
  // ( m1 m2 m3 )
      
	double j1 = A[0];
	double j2 = A[1]; 
	double J = A[2];
	double m1 = A[3];
	double m2 = A[4];
	double M = A[5]; 
    
    double B[6];

    B[0] = j1; 
    B[1] = j2; 
    B[2] = j1 + j2; //j3
    B[3] = m1; 
    B[4] = m2; 
    B[5] = -(m1 + m2); //m3

  return pow(-1, j1 - j2 - B[5])/sqrt(2*B[2]+1) * CG2(B);

}
double SixJSymbol(double J1, double J2, double J3, double J4, double J5, double J6){

  // coupling of j1, j2, j3 to J-J1
  // J1 = j1
  // J2 = j2
  // J3 = j12 = j1 + j2
  // J4 = j3
  // J5 = J = j1 + j2 + j3
  // J6 = j23 = j2 + j3

  //check triangle condition
  if( J3 < abs(J1 - J2 ) || J1 + J2 < J3 ) return 0;
  if( J6 < abs(J2 - J4 ) || J2 + J4 < J6 ) return 0;
  if( J5 < abs(J1 - J6 ) || J1 + J6 < J5 ) return 0;
  if( J5 < abs(J3 - J4 ) || J3 + J4 < J5 ) return 0;

  double sixJ = 0;

  for( float m1 = -J1; m1 <= J1 ; m1 = m1 + 1){
    for( float m2 = -J2; m2 <= J2 ; m2 = m2 + 1){
      for( float m3 = -J3; m3 <= J3 ; m3 = m3 + 1){
        for( float m4 = -J4; m4 <= J4 ; m4 = m4 + 1){
          for( float m5 = -J5; m5 <= J5 ; m5 = m5 + 1){
            for( float m6 = -J6; m6 <= J6 ; m6 = m6 + 1){

              double f = (J1 - m1) + (J2 - m2) + (J3 - m3) + (J4 - m4) + (J5 - m5) + (J6 - m6);

              double a1 = ThreeJSymbol( J1, -m1, J2, -m2, J3, -m3); // J3 = j12
              double a2 = ThreeJSymbol( J1, m1, J5, -m5, J6, m6); // J5 = j1 + (J6 = j23)
              double a3 = ThreeJSymbol( J4, m4, J2, m2, J6, -m6); // J6 = j23
              double a4 = ThreeJSymbol( J4, -m4, J5, m5, J3, m3); // J5 = j3 + j12

              double a = a1 * a2 * a3 * a4;
              //if( a != 0 ) printf("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f | %f \n", m1, m2, m3, m4, m5, m6, a);

              sixJ += pow(-1, f) * a1 * a2 * a3 * a4;

            }
          }
        }
      }
    }
  }

  return sixJ;
}

 ///        (double j1, double j1, double j3, double Lp, double L, double j2)
 //         j3 = j1 + j2, L = |j1 - j2|, Lp = L + 1  
//double SixJ2(double J1, double J2, double J3, double J4, double J5, double J6){
double SixJ2(double C[6]){

  double J1 = C[0]; // j1
  double J2 = C[1]; // j1 
  double J3 = C[2]; // K  
  double J4 = C[3]; // Lp = L + 1
  double J5 = C[4]; // L = |j1 - j2|
  double J6 = C[5]; // j2

  // coupling of j1, j2, j3 to J-J1
  // J1 = j1
  // J2 = j2
  // J3 = j12 = j1 + j2
  // J4 = j3
  // J5 = J = j1 + j2 + j3
  // J6 = j23 = j2 + j3

  //check triangle condition
  if( J3 < abs(J1 - J2 ) || J1 + J2 < J3 ) return 0;
  if( J6 < abs(J2 - J4 ) || J2 + J4 < J6 ) return 0;
  if( J5 < abs(J1 - J6 ) || J1 + J6 < J5 ) return 0;
  if( J5 < abs(J3 - J4 ) || J3 + J4 < J5 ) return 0;

  double sixJ = 0;

  for( float m1 = -J1; m1 <= J1 ; m1 = m1 + 1){
    for( float m2 = -J2; m2 <= J2 ; m2 = m2 + 1){
      for( float m3 = -J3; m3 <= J3 ; m3 = m3 + 1){
        for( float m4 = -J4; m4 <= J4 ; m4 = m4 + 1){
          for( float m5 = -J5; m5 <= J5 ; m5 = m5 + 1){
            for( float m6 = -J6; m6 <= J6 ; m6 = m6 + 1){

              double C1[6] = {J2,J3,J1,-m2,-m3,-m1};
              double C2[6] = {J5,J6,J1,-m5,m6,m1};
              double C3[6] = {J2,J6,J4,m2,-m6,m4};
              double C4[6] = {J5,J3,J4,m5,m3,-m4};


              double f = (J1 - m1) + (J2 - m2) + (J3 - m3) + (J4 - m4) + (J5 - m5) + (J6 - m6);

              double a1 = TJ2(C1); // J3 = j12
              double a2 = TJ2(C2); // J5 = j1 + (J6 = j23)
              double a3 = TJ2(C3); // J6 = j23
              double a4 = TJ2(C4); // J5 = j3 + j12

              double a = a1 * a2 * a3 * a4;
              //if( a != 0 ) printf("%4.1f %4.1f %4.1f %4.1f %4.1f %4.1f | %f \n", m1, m2, m3, m4, m5, m6, a);

              sixJ += pow(-1, f) * a1 * a2 * a3 * a4;

            }
          }
        }
      }
    }
  }

  return sixJ;
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

