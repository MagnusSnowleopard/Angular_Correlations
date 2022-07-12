#include "global.h"
#include "GUI.h"

int main ( int argc, char** argv){

	HistoGUI gui; 

	double A0E = 134.327;
	double A2E = -11.7874;
	double A4E = 0.760896;

	double step = 0.0001;

	vector<double> Theta;
    vector<double> AD_I;
    double ad_start = 0.;//-3.14159/2;   
    int adpoints = (3.14159)/(step);
    double theta,Iad = 0.;
    double aaa,aab,aac = 0.;
    for(int i = 0; i < adpoints; i++){

        theta = i*step + ad_start;
		
        aaa = A0E;
        aab  = A2E*(1.5 * pow(theta,2) - 0.5);
        aac  = A4E*(35./8. * pow(theta,4) - 30./8. * pow(theta,2)  +  3./8. );

        Iad = aaa + aab + aac;
		
		cout << theta << "\n";
		
        AD_I.push_back(Iad);
        Theta.push_back(theta);
		

    }

	vector<double> dangle; 
	//= {0.785398,1.5708,2.35619};
	dangle.push_back(0.785398);
	dangle.push_back(1.5708);
	dangle.push_back(2.35619);
	dangle.push_back(1.5619);
	vector<double> dydata;
	//= {129.,110.,129.};
	dydata.push_back(129.);
	dydata.push_back(110.);
	dydata.push_back(129.);
	dydata.push_back(115.);
	vector<double> deydata;
	//= {10.,10.,10.};
	deydata.push_back(10.);
	deydata.push_back(10.);
	deydata.push_back(10.);
	deydata.push_back(10.);

	

	gui.SetData(dangle,dydata);
	gui.SetErrors(deydata);
	gui.SetFit(A0E,A2E,A4E);

	gui.Init();
	gui.Loop();
	gui.Close();







	return 0; 
}
