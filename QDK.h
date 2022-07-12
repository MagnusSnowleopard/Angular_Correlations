#include "global.h"

using namespace std;


double QK2(double Energy, double radius, double distance, double thickness){
	
	
	double Qkn = 0.; 

	double E_mev = Energy/1000;
	
	double E_log = log(E_mev);

	double EL1 = E_log;
	double EL2 = pow(E_log,2);
	double EL3 = EL1*EL2;
	double EL4 = pow(EL2,2);
	double EL5 = EL4*EL1;

	double TT = -1.1907 -0.5372*EL1 - 0.0438*EL2 + 0.0218*EL3 + 0.0765*EL4 + 0.0095*EL5;

	double Tau = exp(TT);

//	printf("TLN = %lf\n",TT);
//	printf("Tau = %lf\n",Tau);

	//calulating attenuation angles
	
	
	double Z1 = radius / (distance + thickness);
	double Z2 = radius / distance; 

	double alpha = atan(Z1);
	double gamma = atan(Z2);
	double beta;


	double BL = 0.;
	double BU = alpha;
	double A = 0.;
	double delx1  = (BU-BL)/1000;

//	printf("alpha = %lf\n",alpha);
//	printf("gamma = %lf\n",gamma);
//	printf("delx1 = %lf\n",delx1);


	double sum1,sum2,sum3 = 0.;
	double sum4,sum5,sum6 = 0.;

	double cosb,sinb,secb,cscb,c2,c4,fac1,fac2,ex1,ex2 = 0.;
 	double term1,term2,term3 = 0.;
	double term4,term5,term6 = 0.;
	int J=0;
	
	int loop_length = 1000;

	for(int i = 0; i<=loop_length; i++){
/*		
		if(i > 0 and i < loop_length){
			J = i%2;	
			//printf("\t\ti = %d\nJ=%d\n",i,J);
			if(J==0){A=2.;
			}else {A=4.;}
			beta = BL+i+delx1;
		}else{A=.1;beta = BL+i+delx1;}
*/
		if(i != 0){
			if(i != loop_length){
				J = i%2;
				if(J==0){
					A = 2.;
				}else{A=4.;beta = BL+i*delx1;}
			}else{A=1.0;beta = BL+i*delx1;}
		}else{ A =1.0;beta = BL+i*delx1;}

	//	printf("Beta = %lf\n",beta);

		cosb = cos(beta);
		sinb = sin(beta);
		secb = 1.0/cosb;
		c2 = pow(cosb,2);		
		c4 = pow(cosb,4);
		fac1 = -1 *Tau *thickness *secb;
		ex1 = exp(fac1);


		term1 = 0.5*(3*c2-1)*(1-ex1)*sinb*A*delx1;
		term2 = 0.125*A*(35*c4-30*c2+3)*(1-ex1)*sinb*delx1;
		term3 = A*(1-ex1)*sinb*delx1;

		sum1 = sum1 +term1;
		sum2 = sum2 +term2;
		sum3 = sum3 +term3;

	}
		

	double ans1 = sum1/3;
	double ans2 = sum2/3;
	double ans3 = sum3/3;	
	
	double LB=alpha;
	double UB=gamma;
	double delx2 = (UB-LB)/1000;

	for(int i = 0; i<=loop_length; i++){
/*		
		if(i > 0 and i < loop_length){
			J2 = i%2;	
			//printf("\t\ti = %d\nJ=%d\n",i,J2);
			if(J2==0){B=2.;
			}else {B=4.;}
			beta2 = LB+i+delx2;
		}else{B=.1;beta2 = LB+i+delx2;}
*/		
		if(i != 0){
			if(i != loop_length){
				J = i%2;
				if(J==0){
					A = 2.;
				}else{A=4.;beta = LB+i*delx2;}
			}else{A=1.0;beta = LB+i*delx2;}
		}else{A =1.0;beta = LB+i*delx2;}

	//	printf("Beta1 = %lf\n",beta);

		cosb = cos(beta);
		sinb = sin(beta);
		secb = 1.0/cosb;
		cscb = 1.0/sinb;
		c2 = pow(cosb,2);		
		c4 = pow(cosb,4);
		fac2 = -1 *Tau *(radius*cscb -distance*secb);
		ex2 = exp(fac2);

		term4 = 0.5*A*(3*c2-1)*(1-ex2)*sinb*delx2;
		term5 = 0.125*A*(35*c4-30*c2+3)*(1-ex2)*sinb*delx2;
		term6 = A*(1-ex2)*sinb*delx2;

		sum4 = sum4 +term4;
		sum5 = sum5 +term5;
		sum6 = sum6 +term6;

	}
	double	ans4=sum4/3;
	double	ans5=sum5/3;
	double	ans6=sum6/3;

/*
	printf("ans1:%lf\n",ans1*100);
	printf("ans2:%lf\n",ans2*100);
	printf("ans3:%lf\n",ans3*100);
	printf("ans4:%lf\n",ans4*100);
	printf("ans5:%lf\n",ans5*100);
	printf("ans6:%lf\n",ans6*100);
*/	
	double QD2  = (ans1+ans4)/(ans3+ans6);
	double QD4  = (ans2+ans5)/(ans3+ans6);
	printf("--------------\n");
	printf("  QD2 = %lf\n",QD2);
	printf("  QD4 = %lf\n",QD4);
	printf("--------------\n");

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
	
	
	double Qkn = 0.; 

	double E_mev = Energy/1000;
	
	double E_log = log(E_mev);

	double EL1 = E_log;
	double EL2 = pow(E_log,2);
	double EL3 = EL1*EL2;
	double EL4 = pow(EL2,2);
	double EL5 = EL4*EL1;

	double TT = -1.1907 -0.5372*EL1 - 0.0438*EL2 + 0.0218*EL3 + 0.0765*EL4 + 0.0095*EL5;

	double Tau = exp(TT);

//	printf("TLN = %lf\n",TT);
//	printf("Tau = %lf\n",Tau);

	//calulating attenuation angles
	
	
	double Z1 = radius / (distance + thickness);
	double Z2 = radius / distance; 

	double alpha = atan(Z1);
	double gamma = atan(Z2);
	double beta;


	double BL = 0.;
	double BU = alpha;
	double A = 0.;
	double delx1  = (BU-BL)/1000;

//	printf("alpha = %lf\n",alpha);
//	printf("gamma = %lf\n",gamma);
//	printf("delx1 = %lf\n",delx1);


	double sum1,sum2,sum3 = 0.;
	double sum4,sum5,sum6 = 0.;

	double cosb,sinb,secb,cscb,c2,c4,fac1,fac2,ex1,ex2 = 0.;
 	double term1,term2,term3 = 0.;
	double term4,term5,term6 = 0.;
	int J=0;
	
	int loop_length = 1000;

	for(int i = 0; i<=loop_length; i++){
/*		
		if(i > 0 and i < loop_length){
			J = i%2;	
			//printf("\t\ti = %d\nJ=%d\n",i,J);
			if(J==0){A=2.;
			}else {A=4.;}
			beta = BL+i+delx1;
		}else{A=.1;beta = BL+i+delx1;}
*/
		if(i != 0){
			if(i != loop_length){
				J = i%2;
				if(J==0){
					A = 2.;
				}else{A=4.;beta = BL+i*delx1;}
			}else{A=1.0;beta = BL+i*delx1;}
		}else{ A =1.0;beta = BL+i*delx1;}

	//	printf("Beta = %lf\n",beta);

		cosb = cos(beta);
		sinb = sin(beta);
		secb = 1.0/cosb;
		c2 = pow(cosb,2);		
		c4 = pow(cosb,4);
		fac1 = -1 *Tau *thickness *secb;
		ex1 = exp(fac1);


		term1 = 0.5*(3*c2-1)*(1-ex1)*sinb*A*delx1;
		term2 = 0.125*A*(35*c4-30*c2+3)*(1-ex1)*sinb*delx1;
		term3 = A*(1-ex1)*sinb*delx1;

		sum1 = sum1 +term1;
		sum2 = sum2 +term2;
		sum3 = sum3 +term3;

	}
		

	double ans1 = sum1/3;
	double ans2 = sum2/3;
	double ans3 = sum3/3;	
	
	double LB=alpha;
	double UB=gamma;
	double delx2 = (UB-LB)/1000;

	for(int i = 0; i<=loop_length; i++){
/*		
		if(i > 0 and i < loop_length){
			J2 = i%2;	
			//printf("\t\ti = %d\nJ=%d\n",i,J2);
			if(J2==0){B=2.;
			}else {B=4.;}
			beta2 = LB+i+delx2;
		}else{B=.1;beta2 = LB+i+delx2;}
*/		
		if(i != 0){
			if(i != loop_length){
				J = i%2;
				if(J==0){
					A = 2.;
				}else{A=4.;beta = LB+i*delx2;}
			}else{A=1.0;beta = LB+i*delx2;}
		}else{A =1.0;beta = LB+i*delx2;}

	//	printf("Beta1 = %lf\n",beta);

		cosb = cos(beta);
		sinb = sin(beta);
		secb = 1.0/cosb;
		cscb = 1.0/sinb;
		c2 = pow(cosb,2);		
		c4 = pow(cosb,4);
		fac2 = -1 *Tau *(radius*cscb -distance*secb);
		ex2 = exp(fac2);

		term4 = 0.5*A*(3*c2-1)*(1-ex2)*sinb*delx2;
		term5 = 0.125*A*(35*c4-30*c2+3)*(1-ex2)*sinb*delx2;
		term6 = A*(1-ex2)*sinb*delx2;

		sum4 = sum4 +term4;
		sum5 = sum5 +term5;
		sum6 = sum6 +term6;

	}
	double	ans4=sum4/3;
	double	ans5=sum5/3;
	double	ans6=sum6/3;

/*
	printf("ans1:%lf\n",ans1*100);
	printf("ans2:%lf\n",ans2*100);
	printf("ans3:%lf\n",ans3*100);
	printf("ans4:%lf\n",ans4*100);
	printf("ans5:%lf\n",ans5*100);
	printf("ans6:%lf\n",ans6*100);
*/	
	double QD2  = (ans1+ans4)/(ans3+ans6);
	double QD4  = (ans2+ans5)/(ans3+ans6);
//	printf("--------------\n");
//	printf("  QD2 = %lf\n",QD2);
//	printf("  QD4 = %lf\n",QD4);
//	printf("--------------\n");

//Now output a file that contains R, D , T , gamma energy, attentuation coeff, q2 and q4
/*
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
*/

	return QD4;
}























