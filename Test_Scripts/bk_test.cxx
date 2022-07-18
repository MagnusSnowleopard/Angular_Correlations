#include "global.h"
using namespace std;
int main(int argc, char** argv){

    double j1 = 2.;
    double j2 = 1.;
    double Sigma = 0.5;

    double j12 = 2*j1;
    double j14 = 4*j1;

    double sigsq = pow(Sigma,2);
    double sum1 = 0.; 

    double am1,amsq,x,ex = 0.;
    //-------------------SOMETHING WRONG WITH CN1.
    for(int i = 0; i <= j14; i = i + 2){
        am1 = 0.5*(i - j12);
        amsq = pow(am1,2);
        x = - (amsq/(2*sigsq));
        ex = exp(x);

        sum1 = sum1 + ex;
 //       cout << "sum 1 = " << sum1 << "\n";
    }
    double cn1 = 1./sum1;	

    cout << "cn1"  << " " << cn1 << "\n"; 
    
    double am11,amsq1,x1,ex1 = 0.;
    //calculate Bk(j1) for gaussian W(m1) or non zero Sigma

    double Bk11,Bk12 = 0.;

    //NEEDS VALUE FIX
    if(Sigma != 0){
        double A[6];
        A[0] = j1; 
        A[1] = j1; 
        A[5] = 0;

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
/*           
                cout << "---------\n";
                cout << "i = " << i << "\n";
                cout << "m = " << m << "\n";
                
                cout << "AM1 = " << am11 << "\n";       
                cout << "AMSQ = " << amsq1 << "\n";       
                cout << "X = " << x1 << "\n";       
                cout << "I = " << II << "\n";       
                
       
         
                cout << "---------\n";
*/         
         
                cgg = 0.5;  

                ex1 = exp(x1);
                tTerm = cn1*ex1*pow(-1,II) * sfact*cgg;
/*        
                cout << "---------\n";
                cout << "i = " << i << "\n";
                cout << "m = " << m << "\n";
                cout << "ex1 = " << ex1 << "\n";
                cout << "tTerm = " << tTerm << "\n";
*/
                cout << "---------\n";
                if(i == 2){ 
                
                Bk11 = Bk11 + tTerm;
                cout << "---------\n";
                cout << "i = " << i << "\n";
                cout << "m = " << m << "\n";
                cout << "tTerm = " << tTerm << "\n";
                cout << "BK1 =" << Bk11 << "\n";
                cout << "---------\n";
                }
                else if(i == 4){ 
                Bk12 = Bk12 + tTerm; 
                
                cout << "---------\n";
                cout << "i = " << i << "\n";
                cout << "m = " << m << "\n";
                cout << "tTerm = " << tTerm << "\n";
                cout << "BK2 =" << Bk12 << "\n";
                cout << "---------\n";


                }
            }
        }
        cout << "Bk11" << " " << Bk11 << "\n";
        cout << "Bk12" << " " << Bk12 << "\n";


    }

    return 0; 

}
