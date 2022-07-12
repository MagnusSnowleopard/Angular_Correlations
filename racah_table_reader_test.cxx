#include "global.h"
#include "Functlib.h"

using namespace std;


int main(int argc,char** argv){
	

	double j1 = 1.5;
	double j2 = 3.5;


	vector<vector<string> > content;
	vector<string> row;
	string line, word; 
	
	vector<string> j1data;
	vector<string> j2data;
	vector<string> Rk20data; 
	vector<string> Rk21data; 
	vector<string> Rk22data; 
	vector<string> Rk40data; 
	vector<string> Rk41data; 
	vector<string> Rk42data; 
		

	fstream file("RkTable.csv");
	if(file.is_open()){
		while(getline(file,line)){
			row.clear();
	
			stringstream str(line);

			while(getline(str, word, ','))
			
				row.push_back(word);
				content.push_back(row);
				j1data.push_back(row[0]);
				j2data.push_back(row[1]);
				Rk20data.push_back(row[2]);			
				Rk21data.push_back(row[3]);			
				Rk22data.push_back(row[4]);			
				Rk40data.push_back(row[5]);			
				Rk41data.push_back(row[6]);			
				Rk42data.push_back(row[7]);			
				 
		}
		cout<< " --- Angular Data File Loaded --- \n";
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
		//			cout << *dit << " " << *dit2 << " Index Equals =  " <<distance(j2datad.begin(),dit2) << "\n";
					nn = distance(j2datad.begin(),dit2);
				}
			}
		}
	}
	
		cout << j1<< "  " << j2 << "  "<< nn << "\n";



/*	for(dit2 = j2datad.begin();dit2 < j2datad.end();dit2++){
		if(*dit2 == j2){
			num2 = dit2; 
		}
		cout << num2<< "\n";
	}
	
*/

	// table is j1 j2 rk20 rk21 rk22 rk40 rk41 rk42
/*
	for(int i = 0; i< 10;i++){
		cout << " "<< j1data[i] << " "<< j2data[i] << " "<< Rk20data[i] << " "<< Rk21data[i] << " "<< Rk22data[i] << " "<< Rk40data[i] << " "<< Rk41data[i] << " "<< Rk42data[i] << "\n";
	}

*/


	return 0; 
}










