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

#include "GUI.h"
#include "Coeff.h"
#include "QDK.h"
#include "Racah.h"
#include "Cleb.h"


//To compile : g++ AD.cxx -o {Input Executable Name} -lX11
//For more information about angular distributions, read Rose and Brinks, and Frank Moore (1988). 

using namespace std;
void menu(){

	std::cout <<"==========================================\n";
	std::cout<< "\t\t Menu Options \t \n";
	std::cout<< "0 - Reading and Instructions \n"; 
	std::cout<< "1 - Plot Chi Sqr \n";
	std::cout<< "2 - Plot Ang.Dis Fit \n"; 
	std::cout<< "3 - Clear J1 J2 Memory \n";
	std::cout<< "4 - Exit \n";
	std::cout<< "5 - Legfit to Data \n"; 
	std::cout<< "6 - Generate Data \n"; 
	std::cout<<"=========================================== \n";} 
void Readme(){
	std::cout<< " The program calculates Chi-Squared values \n";
	std::cout<< " from experimental angular distributions as \n";
	std::cout<< " a function of multipole ratios using the \n"; 
	std::cout<< " theoretical angular distribution formulae \n"; 
	std::cout<< " in rose and brink \n";
	std::cout<< " \n"; 
	std::cout<< " Follow the prompt in order to correctly display \n";
	std::cout<< " \n";
	std::cout<< " To close the gui, press any button. \n";
	std::cout<< " To zoom in, left click then drag and let go. To unzoom\n";
	std::cout<< " press the space bar. To draw, right click\n";
	std::cout<< "                                          \n";
	std::cout<< " ad.txt is generated with geometric stats.\n";
	std::cout<< " Then at the end it generates CG, W3J & W6J\n";
	
}
#ifndef __CINT__ 

int main(int argc,char **argv){
	HistoGUI gui;
	//	menu();
	double j1, j2;


	//detector data input
	double detradius, targetdistance, detthickness = -1;
	double Energy;
	//input .csv file
	double Sigma; 	
	double Feeding; 

//--------------------------------

	detradius = 3.;
	targetdistance = 4.;
	detthickness = 5.;
	
	Energy = 1062.;
	Sigma = .1; 
	Feeding = .1;

	j1 = 2.; 
	j2 = 1.;
//--------------------------------
	
//TEST MODE? y : 1 || n = 0

	int test = 1; 

	int det_param_token = -1;
	int gamma_energy_token = -1;
	int ang_file_token = -1; 
	int sigma_token = -1; 
	int feeding_token = -1;
	int j1j2token = -1; 
 if(test == 1){

	 det_param_token = 1;
	 gamma_energy_token = 1;
	 ang_file_token = -1; 
	 sigma_token = 1; 
	 feeding_token = 1;
	 j1j2token = 1; 

}
/*
	printf("Input detector radius, target distance, and detector thickness : ");
	scanf("%lf,%lf,%lf",&detradius,&targetdistance,&detthickness);


	if(detradius > 0 && targetdistance > 0 && detthickness > 0){
		det_param_token = 1;
	}else{

		do{
			if(detradius< 0|| targetdistance < 0|| detthickness < 0){
				printf("Negative numbers are not allowed!\nRe-enter radius, distance, and thickness : ");
				scanf("%lf,%lf,%lf", &detradius,&targetdistance,&detthickness);
			}
		}while(detradius< 0|| targetdistance < 0|| detthickness < 0);}

	if(detradius > 0 && targetdistance > 0 && detthickness > 0){
		det_param_token = 1; printf(" --- Detector characteristics loaded --- \n");
	}

	printf("Enter the gamma-ray energy : ");
	scanf("%lf", &Energy);
	if(Energy > 0){
		gamma_energy_token = 1;
	}else{
		do{
			if(Energy <= 0){		
				printf("Negative numbers are not allowed!\nRe-enter gamma-ray energy : ");
			scanf("%lf", &Energy);}

		}while(Energy <= 0);}

	if(Energy > 0){
		gamma_energy_token = 1; printf(" --- Gamma Energy loaded --- \n");
	}*/

	QKN(Energy, detradius, targetdistance, detthickness);
	//then calucate QD2 and QD4, and replace the 0 with them. 
	double QD2 = .7;
	double QD4 = .3;

	//input data file of Angular Data (theta, Yexp, Yerr)
	string fname;
	cout<<"Enter the file name : ";
	cin>>fname;

	vector<vector<string> > content;
	vector<string> row;
	string line, word; 
	
	vector<string> adata;
	vector<string> ydata;
	vector<string> eydata; 
	int num = 0;
	fstream file(fname.c_str());
	if(file.is_open()){
		while(getline(file,line)){
			row.clear();
	
			stringstream str(line);

			while(getline(str, word, ','))
			
				row.push_back(word);
				content.push_back(row);
				adata.push_back(row[0]);
				ydata.push_back(row[1]);
				eydata.push_back(row[2]);			
				 
		}
		cout<< " --- Angular Data File Loaded --- \n";
	}else cout<< "Could not open the file\n";
	
	for(int i = 0; i< adata.size();i++){
//		cout << adata[i] << "\n";
//		cout << ydata[i] << "\n";
//		cout << eydata[i] << "\n";
	}
	
	// Legendre Polinomial fit using adata, ydata, eydata;

	vector<int> int_adata(adata.size());
std::transform(adata.begin()+1, adata.end(), std::back_inserter(int_adata), StrToInt);

	vector<double> double_angle(int_adata.size() -4);
	for(int i = 4; i<  int_adata.size(); i++){
	//	printf("angle = %d\n",int_adata[i]);
		double_angle.push_back((double)int_adata[i]*3.14159/180);	
	}

	vector<int> int_ydata(ydata.size());
std::transform(ydata.begin()+1, ydata.end(), std::back_inserter(int_ydata), StrToInt);

	for(int i = 4; i<  int_ydata.size(); i++){
	//	printf("Ydata = %d\n",int_ydata[i]);
	}

	vector<double> double_ydata(int_ydata.size() -4);
	for(int i = 4; i<  int_ydata.size(); i++){
	//	printf("angle = %d\n",int_adata[i]);
		double_ydata.push_back((double)int_ydata[i]);	

	}
	
	vector<int> int_eydata(eydata.size());
std::transform(eydata.begin()+1, eydata.end(), std::back_inserter(int_eydata), StrToInt);

	for(int i = 4; i<  int_eydata.size(); i++){
	//	printf("Error data = %d\n",int_eydata[i]);
	}

	
	for(int i = 3; i< double_angle.size();i++){
	//	printf("rad angle = %lf\n",double_angle[i]);
	
		//displays the angles inputed as doubles in radians. 
	}
	
	for(int i = 3; i< double_ydata.size();i++){
	//	printf("y_data = %lf\n",double_ydata[i]);
	
		//displays the angles inputed as doubles in radians. 
	}
	/*
	for(int i=0; i<content.size(); i++){
		for(int j=0; j<content[i].size(); j++){
			cout<< content[i][j]<<" ";
		}
		cout<<"\n";
	}
	*/
	//legendre fitting 
	// F[xi] = A0 + A2(P2(xi)) + A4(P4(xi))
	//xi data is angle data in radians. 
	//yi data is y-intensity data. 
	
	//to change the number of angles, just mess with the indexing of the read in and arrays. 


	//12 constants. this will build our matrix for gaussian elinmination. 
	//[a1 a2   a3] [A0] = [a4]
	//[a5 a6   a7] [A2] = [a8]
	//[a9 a10 a11] [A4] = [a12]
	
	double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12 = 0.;
	
	for(int i = 3; i< double_angle.size(); i++){
		
//eq1
		a1 += 1; //how many data points, not intensional but convientient.  
		a2 += (1.5 * pow(double_angle[i],2) - 1);
		a3 += (35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. );
		a4 += double_ydata[i]; 
//eq2 
		a5 += (1.5 * pow(double_angle[i],2) - 1);
		a6 += (1.5 * pow(double_angle[i],2) - 1)*(1.5 * pow(double_angle[i],2) - 1);
		a7 += (1.5 * pow(double_angle[i],2) - 1)*(35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. );
		a8 += (1.5 * pow(double_angle[i],2) - 1)*double_ydata[i];
//eq4
		a9 += (35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. );
		a10 += (1.5 * pow(double_angle[i],2) - 1)*(35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. );
		a11 +=(35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. )* (35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. );
		a12 += (35./8. * pow(double_angle[i],4) - 30./8. * pow(double_angle[i],2)  +  3./8. )*double_ydata[i];

		 
	}
/*
	printf("a1 = %lf\n",a1);	
	printf("a2 = %lf\n",a2);	
	printf("a3 = %lf\n",a3);	
	printf("a4 = %lf\n",a4);	
	printf("a5 = %lf\n",a5);	
	printf("a6 = %lf\n",a6);	
	printf("a7 = %lf\n",a7);	
	printf("a8 = %lf\n",a8);	
	printf("a9 = %lf\n",a9);	
	printf("a10 = %lf\n",a10);	
	printf("a11 = %lf\n",a11);	
	printf("a12 = %lf\n",a12);	
*/
//here number of equations are 3 because its a 3 coeff. fit. 
	int n = 3; 
	double mat[3][4] = {{a1,a2,a3,a4},{a5,a6,a7,a8},{a9,a10,a11,a12}};
//where A0,A2,A4 are stored
	double residual[n];

	for(int i = 0; i<n;i++){
		for( int j = i+1; j<n; j++){
			if(abs(mat[i][i]) < abs(mat[i][j])){
				for(int k = 0; k<n+1;k++){
					//swapping matrix
		//			mat[i][k] = mat[i][k] + mat[j][k];
		//			mat[j][k] = mat[i][k] - mat[j][k];
		//			mat[i][k] = mat[i][k] - mat[j][k];
					swap(mat[i][k],mat[j][k]);
				}

			}
		}
	}
//performing gaussian elimination.

	for(int i = 0; i< n-1; i++) {

		for(int j = i +1; j<n;j++){
			
			float f = mat[j][i]/mat[i][i];

			for(int k = 0; k<n+1; k++){
				mat[j][k] = mat[j][k] -f*mat[i][k];

			}
		}
	}
	
//backwards substitution for unknowns.

	for( int i = n-1; i >= 0; i--){

		residual[i] = mat[i][n];

		for( int j = i + 1; j<n;j++){

			if(i != j){
				residual[i] = residual[i] - mat[i][j] * residual[j];
			}
		}
		residual[i] = residual[i]/mat[i][i]; 
	}
	

	cout<<"A0 = "<<residual[0]<<"\n";
	cout<<"A2 = "<<residual[1]<<"\n";
	cout<<"A4 = "<<residual[2]<<"\n";

	


	/*

	printf("Please enter sigma, sigma = 0 for perfect allignment : ");
	scanf("%lf", &Sigma);
	if(Sigma >= 0 ){
		sigma_token = 1; 
	}else{
		do{

			printf("Negative numbers are not allowed!\nRe-enter Sigma : ");
			scanf("%lf", &Sigma);
		}while(Sigma <0);}

	if(Sigma >= 0 ){
		sigma_token = 1; printf(" --- Sigma Loaded --- \n"); 
	}

	printf("Please enter feeding parameter : ");
	scanf("%lf", &Feeding);	
	if(Feeding >= 0 and Feeding < 1){
		feeding_token = 1;
	}else{
		do{
			printf("Negative numbers are not allowed!\nRe-enter Feeding : ");
			scanf("%lf", &Feeding);
		}while(Feeding < 0);}

	if(Feeding >= 0 ){
		feeding_token = 1; printf(" --- Feeding Loaded --- \n");
	}

	printf("Please enter J1,J2 : ");
	scanf("%lf,%lf",&j1,&j2);	
	if(j1 >= 0 && j2  >= 0 ){
		j1j2token = 1;
	}else{
		do{
			printf("Negative numbers are not allowed!\nRe-enter J1,J2 : ");
			scanf("%lf,%lf", &j1,&j2);	
		}while(j1 < 0 || j2 < 0);}

	if(j1 >= 0 && j2  >= 0 ){
		j1j2token = 1; printf(" --- J1 & J2 Loaded --- \n");
	}
*/
	

	/*
	   printf("detector radius = %lf\n",detradius);
	   printf("detector distance = %lf\n",targetdistance);
	   printf("detector width = %lf\n",detthickness);

	   printf("Gamma-Ray energy = %lf\n", Energy);
	   printf("Sigma Distribution = %lf\n", Sigma); 
	   printf("Feeding param = %lf\n", Feeding); 
	   printf("J1,J2 = %lf,%lf\n",j1,j2);
	   */
/*	
	double J  = j1 + j2;
	double m1 = j1;
	double m2 = j2;
	double M  = m1 + m2; 
	double j3 = J; 
	double m3 = j3;	
	double j4 = j3;
	double j5 = j1 + j2 + j3;
	double j6 = j2 + j3; 
*/
// Calculate BK(j1) for perfect allignment
//
	double A[6];
	double js = sqrt(2*j1 +1); 
	double j12 = 2*j1;
	int IS = (int)j12 % 2; 
	
	double cg12,cg22,cg24,cg14 = 0.;
	double I1, I2 = 0.;
	double Bk11, Bk12 = 0.;
	
	double cg1, cg2 = 0.;	
	
	if(IS == 1){
		A[0] = j1;
		A[1] = j1;
		A[2] = 2.;
		A[3] = .5;
		A[4] = -.5;
		A[5] = 0.;
	
		cg12 = CG2(A);
		
		A[3] = -.5;
		A[4] = .5;
		
		cg22 = CG2(A);

		A[2] = 4.;
		
		cg24 = CG2(A);

		A[3] = .5;
		A[4] = -.5;

		cg14 = CG2(A);

		Bk11 = (0.5)*js*(pow(-1,I1) * cg12 + pow(-1,I2) * cg22);
		Bk12 = (0.5)*js*(pow(-1,I1) * cg14 + pow(-1,I2) * cg24);
	}else if(IS == 0){
		A[0] = j1; 
		A[1] = j1; 
		A[2] = 2.;
		A[3] = 0.;
		A[4] = 0.;
		A[5] = 0.; 
		
		cg1 = CG2(A);
	
		Bk11 = pow(-1,j1) * js * cg1;

		A[2] = 4.;
	
		cg2 = CG2(A);
		Bk12 = pow(-1,j1) * js * cg2; 	
		
	}

//normalize W(m1) 
//

	       j12 = 2*j1;
	double j14 = 4*j1;
	
	double sigsq = pow(Sigma,2);
	double sum1 = 0.; 

	double am1,amsq,x,ex = 0.;

	for(int i = 0; i <= j14; i = i + 2){
		am1 = 0.5*(i - j12);
		amsq = pow(am1,2);
		x = - (amsq/(2*sigsq));
		ex = exp(x);

		sum1 = sum1 + ex;
	}
	double cn1 = 1./sum1;	

	double AL0 = 0.;
	double A0 = j1 - j2; 
	double A0b = abs(A0);
	if(A0b > 0){ AL0 = A0b;}
	
	if(A0b <= 0){AL0 = 1.;}

	double AL1 = AL0 +1; 

	double am11,amsq1,x1,ex1 = 0.;
//calculate Bk(j1) for gaussian W(m1) or non zero Sigma

	A[0] = j1; 
	A[1] = j1; 
	
	double sfact = sqrt(2*j1 +1); 
	double cgg = 0.; 
	double tTerm = 0.; 
	double II = 0.;
	for(int i = 2; i<= 4; i = i + 2){
		
		A[2] = i; 
		Bk11 = 0.;
		Bk12 = 0.; 
		for(int m = 0; m <= j14; m = m +2){
			am11 = 0.5*(m - j12);
			amsq1 = pow(am11,2);
			x1 = -(amsq/(2*sigsq));
			II = j1 - am11; 
			A[3] = am11; 
			A[4] = -am11; 
			cgg = CG2(A);
			ex1 = exp(x1);
			tTerm = cn1*ex1*pow(-1,II) * sfact*cgg;
			if(i == 2) Bk11 = Bk11 * tTerm;
			if(i == 4) Bk12 = Bk12 * tTerm; 
		}
	}
	

//need to add correction of statistical tensors for cascade feeding. 
//
//


//Calculate Rk(LL,j1j2)
//
	double B[6];
	
	double c0,c1,c2,q0,q1,q2,i0,i1,i2,sf0,sf1,sf2, rk01,rk11,rk21, rk02, rk12,rk22 = 0.; 
	
	for(int k = 1; k < 4; k = k + 2){

		A[0] = AL0;
		A[1] = AL0;
		A[2] = k; 
		A[3] = 1.;
		A[4] = -1.;
		A[5] = 0.;
		
		B[0] = j1; 
		B[1] = j1; 
		B[2] = AL0; 
		B[3] = AL0; 
		B[4] = k; 
		B[5] = j2; 
	
		c0 = CG2(A);
		q0 = RACAH(B);
		i0 = 1+j1 -j2-k; 

		sf0 = sqrt((2*j1+2) * (2*AL0 +1) * (2*AL1 +1)); 

		if(k == 2){ rk01 = (pow(-1,i0) * sf0 * c0 * q0);}	
		if(k == 4){ rk02 = (pow(-1,i0) * sf0 * c0 * q0);}

		A[1] = AL1; 
		B[3] = AL1; 

		c1 = CG2(A);
		q1 = RACAH(B);
		i1 = 1 +j1 - j2 + AL1 - AL0 - k;
		sf1 = sqrt((2*j1 +1) * (2*AL0 +1) *(2*AL1 +1));
		 
		if(k == 2){ rk11 = (pow(-1,i1) * sf1 * c1 * q1);}
		if(k == 4){ rk12 = (pow(-1,i1) * sf1 * c1 * q1);}

		A[0] = AL1; 
		B[2] = AL1; 
		
		c2 = CG2(A); 
		q2 = RACAH(B); 
		i2 = 1 + j1 - j2 - k;
		sf2 = sqrt((2*j1 +1) * (2*AL1+1) *(2*AL1+1));

		if(k == 2){ rk21 = (pow(-1,i2) * sf2 * c2 * q2);}
		if(k == 4){ rk22 = (pow(-1,i2) * sf2 * c2 * q2);}
		

	}





	
	printf("Bk1 = %lf\n",Bk11);
	printf("Bk2 = %lf\n",Bk12);

	printf("rk01 = %lf\n",rk01); 
	printf("rk11 = %lf\n",rk11); 
	printf("rk21 = %lf\n",rk21); 
	printf("rk02 = %lf\n",rk02); 
	printf("rk12 = %lf\n",rk12); 
	printf("rk22 = %lf\n",rk22); 
	

/*

	double JJ[6] = {j1,j1,2,.5,-.5,0};
	double racah = RACAH(JJ);	
	double cg222 = CG2(JJ);

	
	
//	double cg = CG(J,M,j1,m1,j2,m2);
//	double W3J= ThreeJsym(j1,m1,j2,m2,j3,m3);
//	double W6J = SixJsym(j1,j2,j3,j4,j5,j6);

	double cg_a = CG_a(JJ);
	double W3J_a= ThreeJsym_a(JJ);
	double W6J_a = SixJsym_a(JJ);
	printf(" CG convert : %lf\n",cg222);	
	printf(" CG         : %lf\n", cg_a);
	printf(" W3-J 	    : %lf\n", W3J_a);
	printf(" W6-J       : %lf\n", W6J_a);	
	printf(" Racah convert  : %lf\n", racah);

*/	
	// here the Chi-sqs from A2 and A4 needed to be added together. 

	 
	double delta_min = -3.14159 / 2.;
	double delta_max = 3.14159 / 2.;
	double step = 0.001; 

 	double delta = 0.;
	double tan_delta =0.;
	double A0E = residual[0];//NEED FrOM AF FIT
	double A2E = residual[1];
	double A4E = residual[2]; 
	double A2T,a2T = 0.;
	double A4T,a4T = 0.;
	double a4E = residual[2]; //NEED FROM AD FIT
	double rd2T = 0.;
	double rd4T = 0.;
	double X22,X24 = 0.;
	double X2_total = 0.;

	int points = (delta_max - delta_min) / step ; 
	
	vector<double> chisqr;
	vector<double> tdelta;

	for(int i = 0; i< points; i++){
		tan_delta = i + step + delta_min;
		delta = tan(tan_delta);
		rd2T = (rk01 + 2*delta*rk11 + pow(delta,2)*rk21)/(1+pow(delta,2));
		A2T = QD2*Bk11*rd2T;
		a2T = A2T/A0E;

		
		rd4T = (rk02 + 2*delta*rk12 + pow(delta,2)*rk22)/(1+pow(delta,2));
		A4T = QD4*Bk12*rd4T;
		a4T = A4T/A0E;
	
		X22 = pow(abs(A2E- a2T),2)/(abs(a2T));
		X24 = pow(abs(A4E -a4T),2)/(abs(a4T));
		X2_total = X22 + X24;

		chisqr.push_back(X2_total);
		tdelta.push_back(tan_delta);

	}






	int param;

	param = param_run(det_param_token, gamma_energy_token, ang_file_token, sigma_token, feeding_token, j1j2token);  	


	int optnum = -1;
	menu();

	printf("Input option : ");
	scanf("%d", &optnum);
	if(optnum < 0){
		optnum = -1;
		printf("Negative numbers are not allowed!\nInput option : ");
		scanf("%d", &optnum);
	}

	if(optnum == 0){
		Readme();
		optnum = -1;

		printf("Input option : ");
		scanf("%d", &optnum);
	} 

	//if all inputs are valid, plot distribution. 
	if(param == 1 && optnum == 1){
		optnum = -1;

		//plotting values here
		std::vector<double> x;
		std::vector<double> y;
		
	        //x will be arctan of mixing ratio. 
	        //y will be log(X^2); 
	        
				

		for(int i = 0; i < 50; i++){
			x.push_back( i );
			y.push_back( pow(i, 0.5) );
		}
		gui.SetData(x,y);
		gui.Init();
		gui.Loop();
		gui.Close();
/*
		printf("Input option : ");
		scanf("%d", &optnum);
		if(optnum < 0){		
			printf("Input option : ");
			scanf("%d", &optnum);
		}*/	
	}



	return 1;


}
#endif














