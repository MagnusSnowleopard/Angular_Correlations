#include <X11/Xlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <cmath> 
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
//To compile : g++ AD.cxx -o {Input Executable Name} -lX11

//#include "Coeff.h"

using namespace std;



double CG2(double A[6]){

	double ii[6];
	double im[6];
	double ix[9];
	double nx[5];

	for(int i = 0;i<6;i++){
		ii[i] = 2.0*A[i];
	}
	double Cleb = 0.;
	double IA = ii[0];
	double IB = ii[1];
	double IC = ii[2];
	double ID = ii[3];
	double IE = ii[4];
	double IF = ii[5];

	if(ID + IE - IF != 0){return 0;}
	double K0 = IA + IB + IC;
	if((int)K0 % 2 !=0){return 0;}

	double K1 = IA + IB - IC;
	double K2 = IC + IA - IB;
	double K3 = IB + IC - IA;
	double K4 = IA - abs(IB-IC);
	double K5 = IB - abs(IC-IA);
	double K6 = IC - abs(IA-IB);

	double karr[6]  = {K1,K2,K3,K4,K5,K6};
	double kmn=karr[0];
	double kmx=karr[0];
	for(int i=0;i<6;i++)
	{
		if(kmn>karr[i])
		{
			kmn=karr[i];
		}
		else if(kmx<karr[i])
		{
			kmx = karr[i];
		}
	}
	double K7 = kmn; 

	printf("K1 = %lf\n",K1);
	printf("K2 = %lf\n",K2);
	printf("K3 = %lf\n",K3);
	printf("K4 = %lf\n",K4);
	printf("K5 = %lf\n",K5);
	printf("K6 = %lf\n",K6);
	printf("K7 = %lf\n",K7);


	if(K7 <  0) return 0; 

	for(int i = 0;i<3;i++){
		if((int)(ii[i]+ii[i+3])%2 !=0)return 0;
		if(ii[i] < abs(ii[i+3]))return 0; 
	}
	double sgnfac =1.0;

	for(int i = 0; i<6;i++){
		im[i] = ii[i];
	}
	double FC2, IT = 0.;
	//-------------------------------------------------------------------------
	if(IA >= IB){
		if(IC >= IB){
			if(IB != 0){
				if(IE < 0){
				//assuming the logic passes;
					FC2 = IC +1;


					IT  = K0/2 + 1;
					ix[0] = IT-IC;
					ix[1] = IT-IB;
					ix[2] = IT-IA;
					ix[3] = (IA+ID)/2 +1;
					ix[4] = ix[3] - ID;
					ix[5] = (IB+IE)/2 + 1;
					ix[6] = ix[5] - IE;
					ix[7] = (IC+IF)/2 +1;
					ix[8] = ix[7] - IF;  
					//lgamma(arg) returns the log of the factorial of arg if arg is a natural number.
					double sqfclg = log(FC2) - lgamma(IT+1);
					double IXI;
					for(int i = 0; i<9;i++){
						IXI = ix[i];
						sqfclg = sqfclg + lgamma(IXI);
					}
					sqfclg = 0.5 * sqfclg;

					double xarr[3] = {ix[0],ix[4],ix[5]};

					double xmn=xarr[0];
					double xmx=xarr[0];
					for(int i=0;i<3;i++)
					{
						if(xmn>xarr[i])
						{
							xmn=xarr[i];
						}
						else if(xmx<xarr[i])
						{
							xmx = xarr[i];
						}
					}
					double NZMAX = xmn;

					double NZ2 = (IB-IC-ID)/2;
					double NZ3 = (IA-IC+IE)/2;

					double carr[3] = {0, NZ2, NZ3};

					double cmn=carr[0];
					double cmx=carr[0];
					for(int i=0;i<3;i++)
					{
						if(cmn>carr[i])
						{
							cmn=carr[i];
						}
						else if(cmx<carr[i])
						{
							cmx = carr[i];
						}
					}


					double NZMIN = cmx+1;

					if(NZMAX > NZMIN){return 0;}
					double SS=0.;
					double S1= pow((-1.),(NZMIN-1));
					for(int NZ = NZMIN;NZ<=NZMAX;NZ++){

						double NZM1 = NZ -1;
						nx[0] = ix[0]-NZM1;
						nx[1] = ix[4]-NZM1;
						nx[2] = ix[5]-NZM1;
						nx[3] = NZ-NZ2;
						nx[4] = NZ-NZ3;
						double termlg = sqfclg - lgamma(NZ);
						for(int i = 0; i<5; i++){
							IXI = nx[i];
							termlg= termlg -lgamma(IXI);
							SS = SS + S1*exp(termlg);
							S1 =-S1;
							Cleb=sgnfac*SS;
						}
					} 
				}else{
	// if (IE >=0) _>>>
					ID = -ID;
					IE = -IE;
					IF = -IF;
					if((int)(IA+IB-IC)/2 %2 != 0) sgnfac = (-1)*sgnfac;

					IT  = K0/2 + 1;
					ix[0] = IT-IC;
					ix[1] = IT-IB;
					ix[2] = IT-IA;
					ix[3] = (IA+ID)/2 +1;
					ix[4] = ix[3] - ID;
					ix[5] = (IB+IE)/2 + 1;
					ix[6] = ix[5] - IE;
					ix[7] = (IC+IF)/2 +1;
					ix[8] = ix[7] - IF;  
					//lgamma(arg) returns the log of the factorial of arg if arg is a natural number.
					double sqfclg = log(FC2) - lgamma(IT+1);
					double IXI;
					for(int i = 0; i<9;i++){
						IXI = ix[i];
						sqfclg = sqfclg + lgamma(IXI);
					}
					sqfclg = 0.5 * sqfclg;

					double xarr[3] = {ix[0],ix[4],ix[5]};

					double xmn=xarr[0];
					double xmx=xarr[0];
					for(int i=0;i<3;i++)
					{
						if(xmn>xarr[i])
						{
							xmn=xarr[i];
						}
						else if(xmx<xarr[i])
						{
							xmx = xarr[i];
						}
					}
					double NZMAX = xmn;

					double NZ2 = (IB-IC-ID)/2;
					double NZ3 = (IA-IC+IE)/2;

					double carr[3] = {0, NZ2, NZ3};

					double cmn=carr[0];
					double cmx=carr[0];
					for(int i=0;i<3;i++)
					{
						if(cmn>carr[i])
						{
							cmn=carr[i];
						}
						else if(cmx<carr[i])
						{
							cmx = carr[i];
						}
					}


					double NZMIN = cmx+1;

					if(NZMAX > NZMIN){return 0;}
					double SS=0.;
					double S1= pow((-1.),(NZMIN-1));
					for(int NZ = NZMIN;NZ<=NZMAX;NZ++){

						double NZM1 = NZ -1;
						nx[0] = ix[0]-NZM1;
						nx[1] = ix[4]-NZM1;
						nx[2] = ix[5]-NZM1;
						nx[3] = NZ-NZ2;
						nx[4] = NZ-NZ3;
						double termlg = sqfclg - lgamma(NZ);
						for(int i = 0; i<5; i++){
							IXI = nx[i];
							termlg= termlg -lgamma(IXI);
							SS = SS + S1*exp(termlg);
							S1 =-S1;
							Cleb=sgnfac*SS;
						}
					} 
					


				}
// if (IB != 0) _>>>>
			}else{
				Cleb = sgnfac;
				return Cleb;

			}
// if (IC < IB) _>>>

		}else{

			IT = IC;	
			IC = IB; 
			IB = IT; 
			IT = IF; 
			IF = -IE;
			IE = -IT;	
			sgnfac = sqrt((2.*A[2]+1.0)/(2.*A[1]+1.0));
			if((int)(im[0]-im[3])/2 %2 != 0){sgnfac = (-1.)*sgnfac;}
			
			Cleb = sgnfac;
			return Cleb; 
		}
	
// if(IA < IB) _>>
	}else{

		if(IA >= IC){
			IT = IC;
			IC = IB;
			IB = IT;
			IT = IF; 	
			IF = -IE; 
			IE = -IT; 			
			
			Cleb = sgnfac; 
			return Cleb; 
// if (IA < IC) as well as (IA < IB)
		}else{
			
			IT =IA;
			IA =IB;
			IB =IT;
			IT =ID;
			ID =IE; 
			IE =IT; 
			
			if((int)(K1/2) %2 != 0) {sgnfac = -1.;}
//call to 135 again.		
			if(IB != 0){
				if(IE < 0){
				//assuming the logic passes;
					FC2 = IC +1;


					IT  = K0/2 + 1;
					ix[0] = IT-IC;
					ix[1] = IT-IB;
					ix[2] = IT-IA;
					ix[3] = (IA+ID)/2 +1;
					ix[4] = ix[3] - ID;
					ix[5] = (IB+IE)/2 + 1;
					ix[6] = ix[5] - IE;
					ix[7] = (IC+IF)/2 +1;
					ix[8] = ix[7] - IF;  
					//lgamma(arg) returns the log of the factorial of arg if arg is a natural number.
					double sqfclg = log(FC2) - lgamma(IT+1);
					double IXI;
					for(int i = 0; i<9;i++){
						IXI = ix[i];
						sqfclg = sqfclg + lgamma(IXI);
					}
					sqfclg = 0.5 * sqfclg;

					double xarr[3] = {ix[0],ix[4],ix[5]};

					double xmn=xarr[0];
					double xmx=xarr[0];
					for(int i=0;i<3;i++)
					{
						if(xmn>xarr[i])
						{
							xmn=xarr[i];
						}
						else if(xmx<xarr[i])
						{
							xmx = xarr[i];
						}
					}
					double NZMAX = xmn;

					double NZ2 = (IB-IC-ID)/2;
					double NZ3 = (IA-IC+IE)/2;

					double carr[3] = {0, NZ2, NZ3};

					double cmn=carr[0];
					double cmx=carr[0];
					for(int i=0;i<3;i++)
					{
						if(cmn>carr[i])
						{
							cmn=carr[i];
						}
						else if(cmx<carr[i])
						{
							cmx = carr[i];
						}
					}


					double NZMIN = cmx+1;

					if(NZMAX > NZMIN){return 0;}
					double SS=0.;
					double S1= pow((-1.),(NZMIN-1));
					for(int NZ = NZMIN;NZ<=NZMAX;NZ++){

						double NZM1 = NZ -1;
						nx[0] = ix[0]-NZM1;
						nx[1] = ix[4]-NZM1;
						nx[2] = ix[5]-NZM1;
						nx[3] = NZ-NZ2;
						nx[4] = NZ-NZ3;
						double termlg = sqfclg - lgamma(NZ);
						for(int i = 0; i<5; i++){
							IXI = nx[i];
							termlg= termlg -lgamma(IXI);
							SS = SS + S1*exp(termlg);
							S1 =-S1;
							Cleb=sgnfac*SS;
						}
					} 
				}else{
	// if (IE >=0) _>>>
					ID = -ID;
					IE = -IE;
					IF = -IF;
					if((int)(IA+IB-IC)/2 %2 != 0) sgnfac = (-1.)*sgnfac;

					IT  = K0/2 + 1;
					ix[0] = IT-IC;
					ix[1] = IT-IB;
					ix[2] = IT-IA;
					ix[3] = (IA+ID)/2 +1;
					ix[4] = ix[3] - ID;
					ix[5] = (IB+IE)/2 + 1;
					ix[6] = ix[5] - IE;
					ix[7] = (IC+IF)/2 +1;
					ix[8] = ix[7] - IF;  
					//lgamma(arg) returns the log of the factorial of arg if arg is a natural number.
					double sqfclg = log(FC2) - lgamma(IT+1);
					double IXI;
					for(int i = 0; i<9;i++){
						IXI = ix[i];
						sqfclg = sqfclg + lgamma(IXI);
					}
					sqfclg = 0.5 * sqfclg;

					double xarr[3] = {ix[0],ix[4],ix[5]};

					double xmn=xarr[0];
					double xmx=xarr[0];
					for(int i=0;i<3;i++)
					{
						if(xmn>xarr[i])
						{
							xmn=xarr[i];
						}
						else if(xmx<xarr[i])
						{
							xmx = xarr[i];
						}
					}
					double NZMAX = xmn;

					double NZ2 = (IB-IC-ID)/2;
					double NZ3 = (IA-IC+IE)/2;

					double carr[3] = {0, NZ2, NZ3};

					double cmn=carr[0];
					double cmx=carr[0];
					for(int i=0;i<3;i++)
					{
						if(cmn>carr[i])
						{
							cmn=carr[i];
						}
						else if(cmx<carr[i])
						{
							cmx = carr[i];
						}
					}


					double NZMIN = cmx+1;

					if(NZMAX > NZMIN){return 0;}
					double SS=0.;
					double S1= pow((-1.),(NZMIN-1));
					for(int NZ = NZMIN;NZ<=NZMAX;NZ++){

						double NZM1 = NZ -1;
						nx[0] = ix[0]-NZM1;
						nx[1] = ix[4]-NZM1;
						nx[2] = ix[5]-NZM1;
						nx[3] = NZ-NZ2;
						nx[4] = NZ-NZ3;
						double termlg = sqfclg - lgamma(NZ);
						for(int i = 0; i<5; i++){
							IXI = nx[i];
							termlg= termlg -lgamma(IXI);
							SS = SS + S1*exp(termlg);
							S1 =-S1;
							Cleb=sgnfac*SS;
						}
					} 
					


				}
// if (IB != 0) _>>>>
			}else{
				Cleb = sgnfac;
				return Cleb;

			}
	







	

		}





	} 



	return Cleb;
}


double CG(double J, double M, double j1, double m1, double j2, double m2){
	//recall that j1,m1 + j2,m2 = J,M

	if(M != m1 + m2) return 0;

	double Jmin  =  abs(j1 - j2);
	double Jmax  =  j1+j2;

	if(J < Jmin || Jmax < J) return 0;
	
	double a0  = (2*J+1.0)*factorial(J+j1-j2) * factorial(J-j1+j2) * factorial(j1+j2-J)/factorial(J+j1+j2+1.0);
	double A0 = sqrt(a0);

	double a = factorial(J+M)   *factorial(J-M);
	double a1= factorial(j1+m1) *factorial(j1-m1);
	double a2= factorial(j2+m2) *factorial(j2-m2);
	double A = sqrt( a  * a1 * a2);

	int pmax = min( min(j1+j2-J,j1-m1),j2 + m2);
	
	double cg = 0.;

	for( int p =0; p<=pmax;p++){

		double p1 = factorial(j1+j2-J-p);
		double p2 = factorial(j1-m1-p);
		double p3 = factorial(j2+m2-p);
		double p4 = factorial(J -j2 +m1 +p);
		double p5 = factorial(J -j1 -m2 +p);
		double t = pow(-1,p)/(factorial(p) * p1 * p2 * p3 * p4 * p5);
			
		cg += t;
	}
	return A0*A*cg;
}

double CG_a(double A[6]){
	//recall that j1,m1 + j2,m2 = J,M
	//testing assingments until the output is the same. 
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
	
	double a0  = (2*J+1.0)*factorial(J+j1-j2) * factorial(J-j1+j2) * factorial(j1+j2-J)/factorial(J+j1+j2+1.0);
	double A0 = sqrt(a0);

	double a = factorial(J+M)   *factorial(J-M);
	double a1= factorial(j1+m1) *factorial(j1-m1);
	double a2= factorial(j2+m2) *factorial(j2-m2);
	double A1 = sqrt( a  * a1 * a2);

	int pmax = min( min(j1+j2-J,j1-m1),j2 + m2);
	
	double cg = 0.;

	for( int p =0; p<=pmax;p++){

		double p1 = factorial(j1+j2-J-p);
		double p2 = factorial(j1-m1-p);
		double p3 = factorial(j2+m2-p);
		double p4 = factorial(J -j2 +m1 +p);
		double p5 = factorial(J -j1 -m2 +p);
		double t = pow(-1,p)/(factorial(p) * p1 * p2 * p3 * p4 * p5);
			
		cg += t;
	}
	return A0*A1*cg;
}
double ThreeJsym(double j1, double m1, double j2, double m2, double j3, double m3){
	//[j1 j2 j3] = (-1)^(j1-j2-m3)/ sqrt(2*j3+1) * CG(j3, -m3, j1, m1, j2, m2)
	//[m1 m2 m3]
	
	
	
	return pow(-1,j1 -j2 -m3)/sqrt(2*j3+1)*CG(j3,-m3,j1,m1,j2,m2);
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

							sixj += pow(-1,h)*b;
	
						}
					}
				}
			}
		}
	}

	return sixj;
}


