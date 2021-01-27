
#ifndef _TRAININGDATA_H
#define _TRAININGDATA_H

#include <iostream>
#include <cstdint>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;


class TrainingData
{
  public:
    
	TrainingData(){ //default constructor
		numSamples = 0;
	}
	
	TrainingData(string fileName, int startcol, int inputsize, int outputsize);
	
	int size(){return numSamples;}
	VectorXi nthInput(int n){
		return inputs.at(n);
	}
	VectorXi nthOutput(int n){
		return outputs.at(n);
	} 
 private:
	int numSamples;
	vector<VectorXi> inputs;
	vector<VectorXi> outputs;
};




#endif
