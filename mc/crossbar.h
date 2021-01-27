
#ifndef _CROSSBAR_H
#define _CROSSBAR_H

#include <iostream>
#include <cstdint>
#include <string>
#include <Eigen/Dense>
#include "spikeneuron.h" 
#include <vector>
#include <time.h>
using namespace std;
using namespace Eigen;


class crossbar
{
  public:
    
	crossbar(){ //default constructor
		vector<int> layers{1};
		crossbar(layers,0,1); 
	}
	crossbar(vector<int> layerBreadth, int reservoirLayer, int _synapseBits); //expects a vector of all of the number of neurons in each layer, and the layer number that the reservoir is on. 0 = no reservoir.
	MatrixXd synapses(){
		return crossbarMatrix;
	}
	MatrixXd mask(){
		return boolmask;
	}
	uint8_t synapseSize(){
		return synapseBits;
	}
	int layerCount(){
		return numLayers;
	}
	int size(){
		return crossbarSize;
	}
	void mutate(int period);
	void operator =(crossbar &rhs){
		crossbarMatrix = rhs.synapses();
		boolmask = rhs.mask();
		synapseBits = rhs.synapseSize();
		score = rhs.score;
		//srand(time(NULL));
	}
	VectorXi actFunct(){
		return actFunctions;
	}
//	uint8_t operator +(uint8_t offset){ 
//		add(offset);
//		return buffer();		
//	}
	int score = 0;    
  private:
	MatrixXd crossbarMatrix;
	MatrixXd boolmask;
	MatrixXi testmat;
	VectorXi actFunctions;
	int synapseBits;
	int crossbarSize;
	int numLayers;
};




#endif
