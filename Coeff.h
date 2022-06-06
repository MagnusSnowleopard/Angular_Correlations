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

using namespace std;


int param_run(int dt, int gt, int at, int st, int ft, int jt){
	int param = 0; 
	if( dt == 1 && gt == 1 && at == 1 && st == 1 && ft == 1 && jt == 1){
		param =0; 
	}else param = 1; 

	return param;
}



int StrToInt(std::string const& s)
{
    std::istringstream iss(s);
    int value;

    if (!(iss >> value)) throw std::runtime_error("invalid int");

    return value;
}


int factorial(int fact){
	for(int i=1;i<=fact;i++){
		fact=fact*i;
	} 
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







