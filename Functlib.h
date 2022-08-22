//################################################################
// This file contains the menu outputs and also misc conversion 
// functions used throughout the program.  
//
//################################################################
#include "global.h"

using namespace std;


void menu(){

	std::cout<< "                                          \n";
	std::cout <<"==========================================\n";
	std::cout<< "\t\t Menu Options \t \n";
	std::cout<< "0 - Reading and Instructions \n"; 
	std::cout<< "1 - Plot Chi Sqr vs arc-tan delta\n";
	std::cout<< "2 - Plot Angular Distribution \n"; 
	std::cout<< "3 - Exit \n";  
	std::cout <<"==========================================\n";
	std::cout<< "                                          \n";
} 

void Readme(){

	std::cout<< "                                            \n";
	std::cout<< " The program calculates Chi-Squared values  \n";
	std::cout<< " from experimental angular distributions as \n";
	std::cout<< " a function of multipole ratios using the   \n"; 
	std::cout<< " theoretical angular distribution formulae  \n"; 
	std::cout<< " in Rose and Brinks.                        \n";
	std::cout<< " \n"; 
	std::cout<< " Follow the prompt in order to correctly display \n";
	std::cout<< " \n";
	std::cout<< " To close the gui, press most buttons.      \n";
	std::cout<< " To zoom in, left click then drag and let go. To unzoom\n";
	std::cout<< " press the space bar. To draw, right click. \n";
	std::cout<< "                                            \n";
	std::cout<< " ad.txt is generated with geometric stats.  \n";
	std::cout<< " Multiple different vectors available for   \n";	
	std::cout<< " other analysis perspectives.               \n";
	std::cout<< "                                            \n";
	std::cout<< " To run in test mode, see the top of AD.cxx.\n";
	std::cout<< " Running in test mode allows for faster use.\n";	
	std::cout<< "                                            \n";
	std::cout<< "                                            \n";
}

int param_run(int dt, int gt, int at, int st, int jt){
	int param = 0; 
	if( dt == 1 && gt == 1 && at == 1 && st == 1  && jt == 1){
		param = 1; 
	}else param = 0; 

	return param;
}



int StrToInt(std::string const& s)
{
    std::istringstream iss(s);
    int value;

    if (!(iss >> value)) throw std::runtime_error("invalid int");

    return value;
}

double StrToDouble(string const& s){
	istringstream iss(s);
	double value;

	if(!(iss >> value)) throw runtime_error("invalid double");

	return value;
}


int factorial(int fact){
	for(int i=1;i<=fact;i++){
		fact=fact*i;
	} 
	return fact;
}

double factorial_d( double fact){

	fact = tgamma(fact + 1);

	return fact;
}


double FACTLOG(int num){

	double faclog[170];
	int RI = 0;

	faclog[0] = 0.;
	faclog[1] = 0.;
	
	for( int i = 3; i < 170; i++){
		RI = i - 1; 
		faclog[i] = log(RI) + faclog[i-1];

	}
	
	double flog = faclog[num];


	return flog;


}


vector<double> string_to_double_vector( vector<string> string_vec){	

	vector<string>::iterator sit;

	vector<double> dv;
		
	if(string_vec.size() < 0 ) throw overflow_error("Invalid String Vector Input\n");
	
	for(sit = string_vec.begin(); sit < string_vec.end(); sit++){

			double d; 

			d = stod(*sit);
	
			dv.push_back(d);

	}
	
	return dv;
}








