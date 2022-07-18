#include "global.h"
//To compile : g++ AD.cxx -o {Input Executable Name} -lX11

//#include "Coeff.h"

using namespace std;

//this program doesnt work and it is not used, I spent so much time transcribing it that maybe, one day some poor soul will come and try to fix it to match the table in rose and brinks. 

double RACAH(double A[6]){
	
	double Racah = 0.;

	double ii[6],ix[6],kk[8];
	
	double nn[4], nx[12],ny[8];

	double K1,K2,K3,K4,K5,K6,K7,K8 = 0.;
	
	double NZMIN, NZMAX = .0;
		
	for(int i = 0;i<6;i++){
		ii[i] = 2.0*A[i];
	}
	
	double IA = ii[0];
	double IB = ii[1];
	double IC = ii[2];
	double ID = ii[3];
	double IE = ii[4];
	double IF = ii[5];

	K1 = IA + IB -IE;
	K3 = IC + ID -IE;
	K5 = IA + IC -IF;
	K7 = IB + ID -IF;
	K2 = IE - abs(IA-IB);
	K4 = IE - abs(IC-ID);
	K6 = IF - abs(IA-IC);
	K8 = IF - abs(IB-ID);


	double karr[8]  = {K1,K2,K3,K4,K5,K6,K7,K8};
	double kmn=karr[0];
	double kmx=karr[0];
	for(int i=0;i<8;i++)
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
	double K0 = kmn; 
	
	if(K0 < 0) return 0;
	
	double nkk, n2kk,nkk2 = 0.;

 	kk[0] = K1;	
 	kk[1] = K2;	
 	kk[2] = K3;	
 	kk[3] = K4;	
	kk[4] = K5;	
 	kk[5] = K6;	
 	kk[6] = K7;	
 	kk[7] = K8;	

	for(int i = 0; i < 7; i = i + 2){
		nkk = kk[i];
		n2kk = kk[i]/2;
		nkk2 = n2kk*2.;
		if(nkk != nkk2){return 0;}

	}
	
	ix[0] = IA;
	ix[1] = ID;
	ix[2] = IB;
	ix[3] = IC;
	ix[4] = IF; 
	ix[5] = IE;
	double Lmin = ix[0];
	int Kmin = 1;
	for( int i = 1; i < 6; i++){
		if(ix[i] - Lmin < 0){
			Lmin=ix[i];
		}
		if(ix[i] - Lmin == 0){
			Kmin=i;
		}
		
		if(ix[i] - Lmin > 0){
			break;
		}
	}
	
	double sgn1 = 1.;
	double KQ,IP,IQ,IT = 0.;
	if(Kmin == 1){
		KQ = 7 - Kmin;
		
		IA = ix[(int)KQ-1];
		ID = ix[Kmin+4-1];
		IE = ix[Kmin-1];
		IF = ix[(int)KQ-4-1]; 

	if(((int)(IA+ID-IE-IF)/2 %2) != 0) {sgn1 = -1.;}

	}else if(Kmin == 2){
	
		KQ = 7 - Kmin;
		
		IA = ix[(int)KQ-1];
		ID = ix[Kmin+4-1];
		IE = ix[Kmin-1];
		IF = ix[(int)KQ-4-1]; 

		if(((int)(IA+ID-IE-IF)/2 %2) != 0) {sgn1 = -1.;}

	}else if(Kmin == 3){
		KQ = 9 - Kmin;
		
		IB = ix[(int)KQ-1];
		IC = ix[Kmin+2-1];
		IE = ix[Kmin-1];
		IF = ix[(int)KQ -2-1];	
		if(((int)(IB+IC-IE-IF)/2)%2 !=0) {sgn1 = -1.;}
	}else if(Kmin == 4){
		
		KQ = 9 - Kmin;
		
		IB = ix[(int)KQ-1];
		IC = ix[Kmin+2-1];
		IE = ix[Kmin-1];
		IF = ix[(int)KQ -2-1];	
		if(((int)(IB+IC-IE-IF)/2)%2 !=0) {sgn1 = -1.;}
	
	}else if(Kmin == 5){

		IE = ix[4];
		IF = ix[5];
		IB = ix[3];
		IC = ix[2];

	}else if(Kmin == 6){
		
	}
	
	IP = IA-IB;
	IQ = IC-ID; 
	
	if(abs(IP) >= abs(IQ)){
		if(IP >= 0){
			IQ = IC-ID;
		}else{

		IT = IB;
		IB = IA; 
		IA = IT; 
		IT = ID; 
		ID = IC; 
		IC = IT; 
		IP = -IP;

		}

	}else{
		IT = IC;
		IC = IA; 
		IA = IT; 
		IT = ID;
		ID = IB;
		IB = IT;
		IP = IQ;
		if(IP >= 0){
			IQ = IC - ID; 
		}else{

		IT = IB;
		IB = IA; 
		IA = IT; 
		IT = ID; 
		ID = IC; 
		IC = IT; 
		IP = -IP;

		}
	}

	double sgn2 = 1.;	
	
	if((int)(IB+ID-IF)/2 %2 != 0) sgn2 = -1.;
	if(IE >= 0){
		double BI = IB; 
		double DI = ID; 
		Racah = sgn1 * sgn2 / sqrt((BI+1.)*(DI+1.));
	}else{
	
	nn[0] = (IA+IB+IE)/2 +1;
	nn[1] = (IC+ID+IE)/2 +1;
	nn[2] = (IA+IC+IF)/2 +1;
	nn[3] = (IB+ID+IF)/2 +1;
	
	nx[0] = nn[0] -IE;
	nx[1] = nn[0] -IB;
	nx[2] = nn[0] -IA;
	nx[3] = nn[1] -IE;
	nx[4] = nn[1] -ID;
	nx[5] = nn[1] -IC;
	nx[6] = nn[2] -IF;
	nx[7] = nn[2] -IC;
	nx[8] = nn[2] -IA;
	nx[9] = nn[3] -IF;
	nx[10] = nn[3] -ID;
	nx[11] = nn[3] -IB;
	
	IP = (IA+IB+IC+ID+4)/2;
	IQ = (IE+IF-IA-ID)/2;
	IT = (IE+IF-IB-IC)/2;

	double a = 1; 
	double b = -IQ+1;
	if(a > b) {

		NZMIN = a; 
	}else{
		NZMIN = b;
	
	}
	
	double kzar[3]  = {nx[1],nx[4],nx[10]};
	double kzmn=kzar[0];
	double kzmx=kzar[0];
	for(int i=0;i<3;i++)
	{
		if(kzmn>kzar[i])
		{
			kzmn=kzar[i];
		}
		else if(kzmx<kzar[i])
		{
			kzmx = kzar[i];
		}
	}
	NZMAX = kzmn;
 
	if(NZMAX < NZMIN){return 0;}
	double sqlog = 0;
	double I1 = 0.;
	for(int i = 0; i <4; i++){
		I1 = nn[i]+1;
		sqlog = sqlog -lgamma(I1);
	}
	for(int i = 0; i <12; i++){
		I1 = nx[i];
		sqlog = sqlog - lgamma(I1);
	}
	sqlog = 0.5*sqlog;

	double sss = 0.;
	double NZM1 = 0.;
	double sslog = 0.;

	for(int nz = NZMIN; nz <= NZMAX; nz++){
		NZM1 = nz - 1; 
		ny[0] = IP - NZM1; 
		ny[1] = nx[0] -NZM1;
		ny[2] = nx[3] -NZM1;
		ny[3] = nx[6] -NZM1;
		ny[4] = nx[9]-NZM1;
		ny[5] = nz;
		ny[6] = IQ+nz;
		ny[7] = IT+nz;

		I1 = ny[0];
		
		sslog = sqlog + lgamma(I1);
		for(int i = 1; i<8; i++){
			I1 = ny[i];
		}sslog = sslog - lgamma(I1);
		sss = sss + (pow(-1.,NZM1)*exp(sslog));
	}
	Racah = sgn1*sss;	

	}
	return Racah;
}
	
	
