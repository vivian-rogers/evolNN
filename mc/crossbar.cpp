// FILE: card_demo.cpp
// This is a small demonstration program showing how the Card and Deck classes are used.
#include <iostream>
#include <cstdint>
//#include <cstdlib>
#include <string>
#include <Eigen/Dense>
#include "spikeneuron.h" 
#include <time.h>
#include <vector>
#include "crossbar.h"
#include <stdlib.h>

using namespace std;
using namespace Eigen;


// PROTOTYPES for functions used by this demonstration program:

crossbar::crossbar(vector<int> layerBreadth, int reservoirLayer, int _synapseBits){ //expects a vector of all of the number of neurons in each layer, and the layer number that the reservoir is on. 0 = no reservoir.
	//goes through and should instantiate eveything
	srand(time(NULL));
	int size = 0;
	for(int i = 0; i < layerBreadth.size(); i++){
		size += layerBreadth[i];
	}
	//crossbarMatrix(size,size);
	//boolmask(size,size);
	crossbarMatrix.setZero(size,size);
	boolmask.setZero(size,size);
	numLayers = layerBreadth.size();
	int row = 0;
	int col = 0;
	cout << "numLayers: " << numLayers << endl;
	for(int i = 0; i < numLayers - 1; i++){ //constructs the boolean mask matrix by accessing submatrices
//		cout << "size: " << size << endl;
		col += layerBreadth[numLayers-i-1];
//		cout << i << "th submatrix, row: col: size: " << row << "  "<< col << "  " << size << endl ;
		boolmask.block(row,col,layerBreadth[numLayers-i-1],layerBreadth[numLayers-i-2]).setOnes(); //fills the submatrix with 1s
		row += layerBreadth[numLayers-i-1];
	}
	synapseBits = _synapseBits;
	crossbarSize = size;
	VectorXi actTemp(size);
	for(int i = 0; i < crossbarSize; i++){
		actTemp(i) = (rand() % 15) - 5;
		if(i < layerBreadth.at(layerBreadth.size() -1)){
			actTemp(i) = 8;	
		}
	}
	actFunctions = actTemp;
	mutate(4);
	//cout << "crossbar size: "<<size <<endl;
}

void crossbar::mutate(int period){
	//time(NULL);
	//cout <<" size " << crossbarSize <<endl;
        for(int i = 0; i < crossbarSize; i++){
		if(rand() % (period) == 0){
			actFunctions(i) += (rand() % 3) - 1;
			if(actFunctions(i) == -8){actFunctions(i)=-7;}
			if(actFunctions(i) == 14){actFunctions(i)=13;}
		}
		for(int j = 0; j < crossbarSize; j++){	
	//		cout << "i: j: " << i << " " <<j <<endl;
			//cout << boolmask << endl <<endl;
	//		cout << crossbarMatrix << endl;
			if((int) boolmask(i,j) == 1){
				if(rand() % period == 0){
					crossbarMatrix(i,j) = rand() % synapseBits;
				}
			}
		}
	}	
}
					
