//#################################################################
//  This program is the main functional of an Angular Ditribution 
//  calculation to fit and determine accepted values of delta.
//
//  Input file need is two .csv, one is a table of racah values, 
//  the second is the experimental data file. The program can handle
//  any amount of angles of detectors.  
//
//  The inputs needed are j1 and j2 values, the gamma-ray  ,the radius, 
//  distance and depth of the detectors, and the sigma of magnetic 
//  substate distribution following the reaction of interest. 
//
//  The program generates 2 graphs - 1 Visible Angular distribution
//  plot and fit normalized by A0. - 2 The log(Chi-squared) vs 
//  arctan(delta) graph. 
//
//  The program also outputs a ad.txt which has QDK coeff. and other
//  detector parameters.  
//
//  To compile : g++ AD.cxx -o {Input Executable Name} -lX11
//  If there is a "stack smashing detection" error add on the flag " -fno-stack-protector "
//For more information about angular distributions, read Rose and Brinks, and Frank Moore (1988). 
//#################################################################

#include "global.h"

#include "GUI_AD.h"
#include "GUI_Base.h"
#include "Functlib.h"
#include "QDK.h"
#include "Cleb.h"


using namespace std;
#ifndef __CINT__ 

int main(int argc,char **argv){
/*
    for( float i = -2.; i < 6.; i += 1.0){
        printf("%.0f! = %f \n",i, factorial(i));
    } 

    for( float m1 = -5.5; m1 <= 5.5; m1 += 1){
        printf("m1: %.1f, %f \n", m1, CGcoeff(2, 0., 5.5, m1, 5.5, -m1));
    }
*/
    double j1, j2;


    //detector data input
    double detradius, targetdistance, detthickness = -1;
    double Energy;
    //input .csv file
    double Sigma; 	 

    //--------------------------------

    detradius = 3.;
    targetdistance = 4.;
    detthickness = 5.;

    Energy = 426.;
    Sigma = 1.0; 


    j1 = 7.5; 
    j2 = 5.5;
    //--------------------------------

    //TEST MODE? y : 1 || n = 0

    int test = 1; 

    int det_param_token = -1;
    int gamma_energy_token = -1;
    int ang_file_token = -1; 
    int sigma_token = -1; 

    int j1j2token = -1; 
    if(test == 1){

        det_param_token = 1;
        gamma_energy_token = 1;
        ang_file_token = -1; 
        sigma_token = 1; 

        j1j2token = 1; 

    }
    if( test == 0 ) {
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
        }

        //then calucate QD2 and QD4, and replace the 0 with them. 
    }

    double QD2 =  QK2(Energy, detradius, targetdistance, detthickness);

    double QD4 =  QK4(Energy, detradius, targetdistance, detthickness);

    //	cout << "QD2 = " << " " << QD2 << "  QD4 = " << " " << QD4 << "\n";

    //input data file of Angular Data (theta, Yexp, Yerr)

    vector<vector<string> > content;
    vector<string> row;
    string line, word; 

    vector<string> adata;
    vector<string> ydata;
    vector<string> eydata; 
    int num = 0;
    string fname,fname2;
    cout<<"Enter the file name : ";
    cin>>fname;

    int aaaa = -1;
    fstream file(fname.c_str());

    if(file.is_open()){
        cout<< " --- Angular Data File Loaded --- \n";
        ang_file_token = 1;
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
    }else{
        do{
            cout<< "Could not open the file\n";
            cout<< "Enter the file name : ";
            cin>>fname2;

            fstream file(fname2.c_str());

            if(file.is_open()){
                aaaa = 1;

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

            }else{ aaaa = -1; }

        }while(aaaa == -1);

        cout<< " --- Angular Data File Loaded --- \n";
        ang_file_token = 1;

    }

    for(int i = 0; i< adata.size();i++){
        //		cout << adata[i] << "\n";
        //		cout << ydata[i] << "\n";
        //		cout << eydata[i] << "\n";
    }

    // Legendre Polinomial fit using adata, ydata, eydata;

    adata.erase(adata.begin());
    ydata.erase(ydata.begin());
    eydata.erase(eydata.begin());

    vector<double> dangle = string_to_double_vector(adata);
    vector<double> dydata = string_to_double_vector(ydata);
    vector<double> deydata = string_to_double_vector(eydata);

    vector<double> dangler;

    for(int i = 0; i < dangle.size(); i++){
        double aa = dangle[i]; 

        dangler.push_back(aa*3.14159/180.);
        printf("Angle %d = %lf\n",i,dangler[i]);
    }
    


    //legendre fitting 
    // F[xi] = A0(1) + A2(P2(cos(xi))) + A4(P4(cos(xi)))
    //xi data is angle data in radians. 
    //yi data is y-intensity data. 

    //12 constants. this will build our matrix for gaussian elinmination. 
    //[a1 a2   a3] [A0] = [a4]
    //[a5 a6   a7] [A2] = [a8]
    //[a9 a10 a11] [A4] = [a12]

    // Note, the y data must be scaled by sin(theta) of the given angle. 

    double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12 = 0.;

    for(int i = 0; i< dangler.size(); i++){

        //eq1
        a1 += 1; //how many data points, not intensional but convientient.  
        a2 += (1.5 * pow(cos(dangler[i]),2) - .5);
        a3 += (35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. );
        a4 += dydata[i]; 
        //eq2 
        a5 += (1.5 * pow(cos(dangler[i]),2) - .5);
        a6 += (1.5 * pow(cos(dangler[i]),2) - .5)*(1.5 * pow(cos(dangler[i]),2) - .5);
        a7 += (1.5 * pow(cos(dangler[i]),2) - .5)*(35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. );
        a8 += (1.5 * pow(cos(dangler[i]),2) - .5)*dydata[i];
        //eq4
        a9 += (35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. );
        a10 += (1.5 * pow(cos(dangler[i]),2) - .5)*(35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. );
        a11 +=(35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. )* (35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. );
        a12 += (35./8. * pow(cos(dangler[i]),4) - 30./8. * pow(cos(dangler[i]),2)  +  3./8. )*dydata[i];


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
    cout<<"a2 = "<<residual[1]/residual[0]<<"\n";
    cout<<"a4 = "<<residual[2]/residual[0]<<"\n";




    if(test == 0){	

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

    }

    // Calculate BK(j1) for perfect allignment
    //

    double A[6];
    double js = sqrt(2*j1 +1); 
    double j12 = 2*j1;
    int IS = (int)j12 % 2; 

    double cg12,cg22,cg24,cg14 = 0.;
    double I1, I2 = 0.;

    double cg1, cg2 = 0.;	

    double Bk10, Bk11, Bk12 = 0.;
    if(Sigma == 0){
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


            printf("cg12 = %lf\n",cg12);
            printf("cg22 = %lf\n",cg22);
            printf("cg24 = %lf\n",cg24);
            printf("cg24 = %lf\n",cg24);

            printf("Bk11 = %lf\n",Bk11);	
            printf("Bk12 = %lf\n",Bk12);	

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

            printf("Bk11 = %lf\n",Bk11);	
            printf("Bk12 = %lf\n",Bk12);	
            printf("cg1 = %lf\n",cg1);
            printf("cg2 = %lf\n",cg2);
        }
    }
    //normalize W(m1) 
    //

    j12 = 2*j1;
    double j14 = 4*j1;

    double sigsq = pow(Sigma,2);
    double sum1 = 0.; 

    double am1,amsq,x,ex = 0.;
    //-------------------------------------------------------SOMETHING WRONG WITH CN1.
    for(int i = 0; i <= j14; i = i + 2){
        am1 = 0.5*(i - j12);
        amsq = pow(am1,2);
        x = - (amsq/(2*sigsq));
        ex = exp(x);

        sum1 = sum1 + ex;
    }
    double cn1 = 1./sum1;	

    //	cout << "cn1"  << " " << cn1 << "\n"; 

    double AL0 = 0.;
    double A0 = j1 - j2; 
    double A0b = abs(A0);
    if(A0b > 0){ AL0 = A0b;}

    if(A0b <= 0){AL0 = 1.;}

    double AL1 = AL0 +1; 

    double am11,amsq1,x1,ex1 = 0.;
    //calculate Bk(j1) for gaussian W(m1) or non zero Sigma

    //NEEDS VALUE FIX
    if(Sigma != 0){
        double A[6];
        A[0] = j1; 
        A[1] = j1; 
        A[5] = 0.;
        double sfact = sqrt(2*j1 +1); 
        double cgg = 0.; 
        double tTerm = 0.; 
        double II = 0.;
        for(int i = 2; i<= 4; i = i + 2){

            A[2] = i; 

            for(int m = 0; m <= j14; m = m +2){
                am11 = 0.5*(m - j12);
                amsq1 = pow(am11,2);
                x1 = -(amsq1/(2*sigsq));
                II = j1 - am11; 
                A[3] = am11; 
                A[4] = -am11; 
                
                

                cgg = CGcoeff(i,0,j1,am11,j1,-am11);
               // cgg = CG2(A);
              //  printf("K : %d, j1 : %.1f m1 : %.1f, %f \n", i, j1, am11, cgg);    
            //    cout << "CG =" << cgg << "\n";
             //   cout << "m1 =" << am11 << "\n";
               // cout << "m2 =" << -am11 << "\n";
                ex1 = exp(x1);
                tTerm = cn1*ex1*pow(-1,II) * sfact*cgg;

                if(i == 2){ Bk11 = Bk11 + tTerm;
                }else if(i == 4){ Bk12 = Bk12 + tTerm;} 
            }
        }

//        cout << "Bk10" << " " << Bk10 << "\n";
        cout << "Bk11" << " " << Bk11 << "\n";
        cout << "Bk12" << " " << Bk12 << "\n";


    }

    vector<vector<string> > content_r;
    vector<string> row_r;
    string line_r, word_r; 

    vector<string> j1data;
    vector<string> j2data;
    vector<string> Rk20data; 
    vector<string> Rk21data; 
    vector<string> Rk22data; 
    vector<string> Rk40data; 
    vector<string> Rk41data; 
    vector<string> Rk42data; 


    fstream file1("RkTable.csv");
    if(file1.is_open()){
        while(getline(file1,line_r)){
            row_r.clear();

            stringstream str(line_r);

            while(getline(str, word_r, ','))

                row_r.push_back(word_r);
            content_r.push_back(row_r);
            j1data.push_back(row_r[0]);
            j2data.push_back(row_r[1]);
            Rk20data.push_back(row_r[2]);			
            Rk21data.push_back(row_r[3]);			
            Rk22data.push_back(row_r[4]);			
            Rk40data.push_back(row_r[5]);			
            Rk41data.push_back(row_r[6]);			
            Rk42data.push_back(row_r[7]);			

        }
        cout<< " --- Rk(llj1j2) Table Loaded --- \n";
    }else cout<< "Could not open the file\n";

    vector<double>::iterator dit;
    vector<double>::iterator dit2;

    vector<double> j1datad   = string_to_double_vector(j1data);
    vector<double> j2datad   = string_to_double_vector(j2data);
    vector<double> rk20datad = string_to_double_vector(Rk20data);
    vector<double> rk21datad = string_to_double_vector(Rk21data);
    vector<double> rk22datad = string_to_double_vector(Rk22data);
    vector<double> rk40datad = string_to_double_vector(Rk40data);
    vector<double> rk41datad = string_to_double_vector(Rk41data);
    vector<double> rk42datad = string_to_double_vector(Rk42data);

    int nn = 0;
    for(dit = j1datad.begin();dit < j1datad.end();dit++){

        if(j1 == *dit){

            for(dit2 = j2datad.begin();dit2 < j2datad.end();dit2++){

                if(j2 == *dit2 and j1 == j1datad[distance(j2datad.begin(),dit2)]){
                    //    			cout << *dit << " " << *dit2 << " Index Equals =  " <<distance(j2datad.begin(),dit2) << "\n";
                    nn = distance(j2datad.begin(),dit2);
                }
            }
        }
    }

    //		cout << j1<< "  " << j2 << "  "<< nn << "\n";


    // rk01,rk11,rk21, rk02, rk12,rk22 //READ OUT OF Rk TABLE 

    double rk01 = rk20datad[nn];
    double rk11 = rk21datad[nn];
    double rk21 = rk22datad[nn];
    double rk02 = rk40datad[nn];
    double rk12 = rk41datad[nn];
    double rk22 = rk42datad[nn];
    //work in progress for possible generation or racah
    /*
       double CC1[6] = {j1,j1,2,abs(j1-j2),abs(j1-j2),j2};
       double CC2[6] = {j1,j1,2,abs(j1-j2),abs(j1-j2)+1,j2};
       double CC3[6] = {j1,j1,2,abs(j1-j2)+1,abs(j1-j2)+1,j2};


       double rkk1 = SixJ2(CC1);
       double rkk2 = SixJ2(CC2);
       double rkk3 = SixJ2(CC3);

       printf("rk1 = %lf\n",rkk1); 
       printf("rk2 = %lf\n",rkk2); 
       printf("rk3 = %lf\n",rkk3); 



       double J1 = C[0]; // j1
       double J2 = C[1]; // j1 
       double J3 = C[2]; // K  
       double J4 = C[3]; // Lp = L + 1
       double J5 = C[4]; // L = |j1 - j2|
       double J6 = C[5]; // j2
       */
    //	printf("Bk1 = %lf\n",Bk11);
    //	printf("Bk2 = %lf\n",Bk12);

    printf("rk01 = %lf\n",rk01); 
    printf("rk11 = %lf\n",rk11); 
    printf("rk21 = %lf\n",rk21); 
    printf("rk02 = %lf\n",rk02); 
    printf("rk12 = %lf\n",rk12); 
    printf("rk22 = %lf\n",rk22); 


    // here the Chi-sqs from A2 and A4 needed to be added together. 


    double delta_min = -3.14159/2.;
    double delta_max =  3.14159/2.;
    double THETA_min = 0.;
    double THETA_max = 3.14159;
    double step = 0.01; 

    double Tangle = 0;
    double delta = 0.;
    double atan_delta =0.;
    double A0E = residual[0];//NEED FrOM AF FIT
    double A2E = residual[1];
    double A4E = residual[2]; 
    double A2T,a2T = 0.;
    double A4T,a4T = 0.;
    double rd0T = 0.;
    double rd2T = 0.;
    double rd4T = 0.;
    int points = (delta_max - delta_min) / step ; 
    double X2_total;
    double YT_tot[points];
    //  vector<double> YTT;
    double a2E = A2E/A0E;
    double a4E = A4E/A0E;


    //int points = (delta_max - delta_min) / step ; 

    vector<double> chisqr;
    vector<double> tdelta;
    vector<double> YT_point;
    vector<double> YE_point;
    double YT0,YT2,YT4,YT,YE = 0.;

    double Y_err; 

    //delta loop;
    for(int i = 0; i< points; i++){
        atan_delta = i*step + delta_min;
        delta = tan(atan_delta);
        //now sum over all theta;

        for(int j = 0; j < dangler.size(); j++){
            Tangle = dangler[j];
            Y_err = deydata[j];

            //calculate Y_intensity_theory and and experiment and sum over theta


            rd2T = (rk01 + 2*delta*rk11 + pow(delta,2)*rk21)/(1+pow(delta,2));

            rd4T = (rk02 + 2*delta*rk12 + pow(delta,2)*rk22)/(1+pow(delta,2));

            YT0 = 1; // k = 0; Bk10*?


            YT2 = QD2*Bk11*rd2T*((1.5*pow(cos(Tangle),2)-.5));

            YT4 = QD4*Bk12*rd4T*(35./8.* pow(cos(Tangle),4) - 30./8.*pow(cos(Tangle),2) + 3./8.);

            YT = YT0 + YT2 + YT4;

            //    YT_tot[i] += YT;

            YE = 1. + a2E*(1.5*pow(cos(Tangle),2)-.5) + a4E*(35./8.* pow(cos(Tangle),4) - 30./8.*pow(cos(Tangle),2) + 3./8.);

            X2_total += pow((YT- YE),2)/((dangler.size()-1)*(pow(Y_err,2)));
            //           X2_total[i] = X2_total[i] + pow((YE),2)/((dangler.size()-1)*(pow(Y_err,2)));
        }
        //   if(i%27 == 1){ 
        //   cout << "X2 = " << X2_total[i] << "\n";
        //   }
        //   YTT.push_back(YT_tot[i]);
        chisqr.push_back(log(X2_total));

        tdelta.push_back(atan_delta);

        X2_total = 0.;

        // Attempt to take 1 delta value and plot YT; 


    }
// this is for viewing the theoretical angular distribution as a function of constant delta.    
       for(int k = 0; k < dangler.size(); k++){
       Tangle = dangler[k];
       delta = .8;

    //calculate Y_intensity_theory and and experiment and sum over theta

    //     rd0T = (1 + 2*delta + pow(delta,2))/(1+pow(delta,2));

    rd2T = (rk01 + 2*delta*rk11 + pow(delta,2)*rk21)/(1+pow(delta,2));

    rd4T = (rk02 + 2*delta*rk12 + pow(delta,2)*rk22)/(1+pow(delta,2));

    YT0 = 1; // k = 0;

    YT2 = QD2*Bk11*rd2T*((1.5*pow(cos(Tangle),2)-.5));

    YT4 = QD4*Bk12*rd4T*(35./8.* pow(cos(Tangle),4) - 30./8.*pow(cos(Tangle),2) + 3./8.);

    YT = YT0 + YT2 + YT4;
    //        cout << "YT = " <<YT << "\n"; 
    YT_point.push_back(YT);

    YE = 1 + a2E*(1.5*pow(cos(Tangle),2)-.5) + a4E*(35./8.* pow(cos(Tangle),4) - 30./8.*pow(cos(Tangle),2) + 3./8.);
    YE_point.push_back(YE);

    //  X2_total += pow((YT- YE),2)/((dangler.size()-1)*(pow(Y_err,2)));
    }
    
    /*
       A2T = QD2*Bk11*rd2T;
       a2T = A2T/A0E;
    //using the idea that A0E is our normalizer. 

    rd4T = (rk02 + 2*delta*rk12 + pow(delta,2)*rk22)/(1+pow(delta,2));
    A4T = QD4*Bk12*rd4T;
    a4T = A4T/A0E;

    X22 = pow((a2E -a2T),2)/(a2T);
    X24 = pow((a4E -a4T),2)/(a4T);
    X2_total = (X22 + X24)/2;
    */
    vector<double> Theta;
    vector<double> AD_I;
    double ad_start = 0.;	
    int adpoints = (3.1415)/(step);
    double theta,Iad = 0.;
    double aaa,aab,aac = 0.;
    for(int i = 0; i < adpoints; i++){

        theta = i*step + ad_start;
        aaa = 1.;
        aab  = (A2E/A0E)*(1.5 * pow(cos(theta),2) - .5);
        aac  = (A4E/A0E)*(35./8. * pow(cos(theta),4) - 30./8. * pow(cos(theta),2)  +  3./8. );

        Iad = aaa + aab + aac;

        AD_I.push_back(Iad);
        Theta.push_back(theta);
        //these arrays arent used, but they are the function (if one just wanted to the legendre fit function that goes over the data. 
    }


    vector<double> dydatas;

    for(int i = 0; i < dydata.size(); i++){
        double bbb = dydata[i]; 

         dydatas.push_back(bbb/A0E);
        //  printf("Normalized y-intensity %i = %lf\n",i+1,dydatas[i]);

    }
    //This is to temporarily fix the issue with the AD distribution, it takes off the last point in the arrays when plotting or doesnt see them. NEED FIX
    dydatas.push_back(A0E);
    dangler.push_back(3.1415/2.);

    int param;

    param = param_run(det_param_token, gamma_energy_token, ang_file_token, sigma_token, j1j2token);  	
    /*	
        cout << det_param_token << "\n"; 
        cout << gamma_energy_token <<  "\n";  
        cout << ang_file_token <<  "\n"; 
        cout << sigma_token <<  "\n";  
        cout << j1j2token <<  "\n"; 

        cout << param << "\n";
        */

    //Initialize Graphics Here. 
    //  cout<< "Dangler size = " << dangler.size() << "\n";
    //  cout<< "YT point size = "<< YT_point.size() << "\n";

    HistoGUI gui;

    HistoGUIad gui_ad; 

    // FIX menu input complexity. 
    //
    int optnum = -1;
    menu();

    printf("Please enter Input Option : ");
    scanf("%d", &optnum);
    if(optnum >= 0){


        // menu
        if(optnum == 0 and param == 1){
            Readme();
            optnum = -1;

            printf("Please enter Input Option : ");
            scanf("%d", &optnum);
            //Chi-squared
        }else if(optnum == 1 and param == 1){

            optnum = -1;
            //plotting values here
            //x will be arctan of mixing ratio. 
            //y will be log(X^2); 

            gui.SetData(tdelta,chisqr);
           // gui.SetData(dangler,YE_point);
           // gui.SetData(dangler,YT_point);
            //gui.SetData(Theta,AD_I);
            gui.Init();
            gui.Loop();
            gui.Close();
            printf("Please enter Input Option : ");
            scanf("%d", &optnum);
        }else if(optnum == 2 and param == 1){

            optnum = -1;
           // gui_ad.SetData(dangler,dydata); //dydata means you  are using already normalizing data.
            gui_ad.SetData(dangler,dydatas); //dydatas means you are using non-A0 normalized data to start with.
            gui_ad.SetErrors(deydata);
            gui_ad.SetFit(residual[0],residual[1],residual[2]);

            gui_ad.Init();
            gui_ad.Loop();
            gui_ad.Close();
            printf("Please enter Input Option : ");
            scanf("%d", &optnum);

        }else if(param != 1){
            cout << "Invalid parameter input, something is wrong\n";

        }


    }else{
        do{
            printf("Negative numbers are not allowed!\nRe-enter Option : ");
            scanf("%d", &optnum);
        }while(optnum < 0);}

    if(optnum >= 0 and optnum != 3){
        //repeat 3 cases.


        // menu
        if(optnum == 0 and param == 1){
            Readme();
            optnum = -1;

            //Chi-squared
        }else if(optnum == 1 and param == 1){

            optnum = -1;
            //plotting values here
            //x will be arctan of mixing ratio. 
            //y will be log(X^2); 

            //gui.SetData(dangler,YT_point);
            gui.SetData(tdelta,chisqr);
            //gui.SetData(Theta,AD_I);
            gui.Init();
            gui.Loop();
            gui.Close();
        }else if(optnum == 2 and param == 1){

            optnum = -1;
            gui_ad.SetData(dangler,dydatas);
            gui_ad.SetErrors(deydata);
            gui_ad.SetFit(residual[0],residual[1],residual[2]);

            gui_ad.Init();
            gui_ad.Loop();
            gui_ad.Close();

        }else if(param != 1){
            cout << "Invalid parameter input, something is wrong\n";

        }


    }



    if(param == 1 and optnum == 3){

        return 0;
    }

    return 1;


}
#endif














