// This is a small demonstration program showing how the Card and Deck classes are used.
#include <iostream>    // Provides cout and cin
#include <cstdlib>     // Provides EXIT_SUCCESS
#include <fstream>
#include <cstdint>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "spikeneuron.h"
#include <time.h> 
#include <vector>
#include "crossbar.h"
#include <stdlib.h>
#include "trainingdata.h"
#include <string.h>
#include <deque>
#include <math.h>
//using Eigen::MatrixXd;
//using Eigen::Dynamic;
using namespace std;
using namespace Eigen;

//typedef spikeNeuron(2000,-500,-50) ourNeuron;
//typedef SparseMatrix<int> SpMatrix;




//typedef Matrix<spikeNeuron, Dynamic, 1> neuronVector;


typedef Matrix<VectorXd, Dynamic, 1> lookupVector;
//typedef Matrix<int, Dynamic, 1> VectorXd_uint8_t;
//typedef Matrix<uint8_t, Dynamic, Dynamic> reservoirXd;


void initPopulation(vector<crossbar> &population, int size);
//void iterate(VectorXd &resV, neuronVector &resN);
void displayPop(vector<crossbar> &population);
VectorXi forwardProp(crossbar &population, VectorXi input);
VectorXi translateInput(crossbar &crossbarArray, VectorXi input);
VectorXi activation(VectorXi resV, VectorXi actVect);
//vector<int> popBest(vector<crossbar> &population, int *score); // returns the best crossbar, returns all of the scores to zero
int costFunction(VectorXi actual, VectorXi expected);
int mutationRate(int generation);
int successfulGuess(VectorXi actual, VectorXi expected);
vector<int> popBest(vector<crossbar> &population, int *score, VectorXi *tempScoresSave); // returns the best crossbar, returns all of the scores to zero

int synapse = 1;
int size;
lookupVector matrixLookup;
const vector<int> layerBreadth{6,19,5,3}; //size of vector in each layer
//const vector<int> layerBreadth{6,8,12,12,6,5,3}; //size of vector in each layer
const int popSize = 50;
const int numSamples = 10;
const int numLeaders = 40;
int mutationSpeed = 4;
VectorXi layerBreadthVect = Map<const VectorXi, Eigen::Unaligned>(layerBreadth.data(), layerBreadth.size());
int main()
{
	srand(time(NULL));
	vector<crossbar> population(popSize,crossbar(layerBreadth,0,2));
//	for(int i = 0; i < popSize; i++){
//		population.at(i).mutate((i % 3) + 1);
//	}
//	MatrixXi testttttt;
	TrainingData data("outfile.csv", 1, layerBreadth.at(0), layerBreadth.at(layerBreadth.size() -1));
	//crossbar test(layerBreadth,0,2);	
//	displayPop(population);	
	//vector<crossbar> population(10, crossbar(layerBreadth, 0, 1));
//	size = population[0].mask().size();	
//	cout << population[0].synapses() << endl;
//	cout << test.mask();	
//	cout << "Setting up " << size << " unit impulse lookup tables... (might take a while)" << endl;
//	std::cout << reservoir << std::endl << std::endl;
//	VectorXd resV(size);
	//neuronVector resN(size);
	//cout << "Setting up " << size << " input tables...";
	//SparseVector<uint8_t> vec(size);
				
//	vector<int> test{0,5,0,0,4,7}; //size of vector in each layer
	bool rpt = true;
	int bestCrossbar;
	int sampleIndex;
	crossbar temp(layerBreadth,0,2);
	VectorXi inputV(layerBreadth.at(0));
	VectorXi outputV(layerBreadth.at(layerBreadth.size() -1));
	int generation = 0;
	int indScore;
	int bestScore;
	int bestIndex;
	int randData;
	int scoreSum = 0;
	int correctSum = 0;
	vector<int> sortedPop;
	VectorXi sortedScores(population.size());
	deque<float> scoreQueue;
	deque<float> correctnessQueue;
	VectorXi numCorrect(population.size());
	int tempNumcorrect;
	while(rpt){
		numCorrect.setZero(population.size());
		cout << endl << "Generation: " << generation << ", Layer sizes: " << layerBreadthVect.transpose() <<", Samples: " << numSamples << ", Mutation period: " << mutationRate(generation) << endl ;
		for(int i = 0; i < numSamples; i++){
	//		sampleIndex = i; 
			sampleIndex = rand() % data.size(); 
		//	cout <<sampleIndex <<",	";
			for(int j = 0; j < popSize; j++){
				//test = {rand()%8,rand()%8,rand()%8,rand()%8,rand()%8,rand()%8};
				inputV = translateInput(population.at(j), data.nthInput(sampleIndex));
	//			cout << population.at(0).synapses() <<endl;
	//			cout << population.at(0).synapses().cast<int>() * input <<endl;
				outputV = forwardProp(population.at(j), inputV);
				indScore = costFunction(outputV, data.nthOutput(sampleIndex));	
				population.at(j).score += indScore;
				numCorrect(j) += successfulGuess(outputV, data.nthOutput(sampleIndex));
	//				cout << population.at(0).mask() <<endl;		
	//			cout << population.at(0).synapses() <<endl;		
	//			cout <<"output: " <<output.transpose() <<endl <<endl;
	//			rpt = false;
			//	for(int i = 0; i <10000000; i++){}
			}
		}
		sortedPop = popBest(population, &bestScore, &sortedScores);
		//temp = population.at(sortedPop.at(0));
		//population.at(0) = temp;
		//population.at(0).score = 0;
		scoreQueue.push_back((float) bestScore/numSamples);
		//scoreQueue.push_back((float) bestScore/(layerBreadth.at(layerBreadth.size()-1)*numSamples));
		if(scoreQueue.size() > 200){scoreQueue.pop_front();}
		if(scoreQueue.size() > 0){
			for(int i = 0; i < scoreQueue.size(); i++){
				scoreSum += scoreQueue.at(i);
			}
		}
		correctnessQueue.push_back((float) (100 * numCorrect(sortedPop.at(0)))/numSamples);
		if(correctnessQueue.size() > 200){correctnessQueue.pop_front();}
		if(correctnessQueue.size() > 0){
			for(int i = 1; i < correctnessQueue.size(); i++){
				correctSum += correctnessQueue.at(i);
			}
		}
		cout << "Moving avg error^2: "<< (float) scoreSum / ( scoreQueue.size() * layerBreadth.at(layerBreadth.size() -1)) << "		";
		//cout << "Moving correct%: "<< (float) correctSum/scoreQueue.size() << "		";
		cout << "Best score: "<< scoreQueue.back() << endl;
		cout << "Sorted Scores:"<< sortedScores.transpose() << endl;
		//cout << "Number correct: {" << numCorrect.transpose() << "}/" << numSamples<<endl <<endl;
		correctSum = 0;
		scoreSum = 0;
		for(int i = numLeaders; i < popSize; i++){
			//if(rand() % 2 == 0){
			population.at(sortedPop.at(i)) = population.at(sortedPop.at(rand() % (numLeaders/3)));
			population.at(sortedPop.at(i)).mutate(mutationRate(generation));
			//}//population.at(i).score = 0;
		}
		if(generation % 100 == 0){
			randData = rand() % data.size();
			cout << "Best crossbar array: " <<endl;
			cout << population.at(sortedPop.at(0)).synapses() << endl <<endl;
			//cout << population.at(sortedPop.at(0)).mask() << endl <<endl;
			cout << "Activation function offsets: " <<population.at(sortedPop.at(0)).actFunct().transpose() << endl;
			cout <<  "Input: " <<data.nthInput(randData).transpose() << " -> ";
			outputV =  forwardProp(population.at(sortedPop.at(0)), translateInput(population.at(sortedPop.at(0)), data.nthInput(randData)));
			cout << outputV.transpose(); 
//			cout << (population.at(sortedPop.at(0)).synapses() * data.nthInput(randData) ).transpose();
			cout << " vs " << data.nthOutput(randData).transpose() <<endl <<endl;
		}
		generation++;
	}
}


int successfulGuess(VectorXi actual, VectorXi expected){
	int highestScoreActual = 0;
	int highestScoreExpected = 0;
	int highestIndexActual = 0;
	int highestIndexExpected = 0;
//	cout << "actual: " <<actual.transpose() << " Expected: "<<expected.transpose() <<endl;
	for(int i=0; i <actual.size(); i++){
		if(actual(i) > highestScoreActual){
			highestScoreActual = actual(i);
			highestIndexActual = i;
		}
		if(expected(i) > highestScoreExpected){
			highestScoreExpected = expected(i);
			highestIndexExpected = i;
		}
	}
	if(highestIndexExpected == highestIndexActual){
		return 1;
	}
	else{
		return 0;
	}
}


int mutationRate(int generation){
	return (sqrt(generation/mutationSpeed + 1000)/2) - 10;	
}

int costFunction(VectorXi actual, VectorXi expected){
	int sum = 0;
	int dif;
	for(int i = 0; i < actual.size(); i++){
		dif = actual(i) - expected(i);
		//if(false){
		//if(successfulGuess(actual,expected)){
		//	if(dif > 0){
		//		sum += -2;
		//		//sum -= dif*dif;
		//	}
		//}else{
		//	if(dif < 0){
		//		sum += (-1) * dif * dif * dif;
		//	}
		//	else{
		//		sum += dif * dif * dif;
		//	}
		sum += dif*dif;
		//} 
		//sum += dif*dif;
	}
	return sum;
}

vector<int> popBest(vector<crossbar> &population, int *score, VectorXi *tempScoresSave){ // returns the best crossbar, returns all of the scores to zero
	vector<int> tempIndices;
	vector<int> tempScores;
	int bestScore = 100000000;
	int bestIndex = 0;
	for(int i = 0; i < popSize; i++){
		tempIndices.push_back(i);
		tempScores.push_back(population.at(i).score);
		population.at(i).score = 0;
		//cout << tempScores.at(i) << "  ";
//		if(population.at(i).score < bestScore && population.at(i).score > 0){
//			bestIndex = i;
//			bestScore = population.at(i).score;
//		}
	//	cout << " " << population.at(i).score << " ";
	}
	//cout << bestScore << "   ";
	int swapScore;
	//cout <<endl << "rank: " ;
	//cout <<"Popsize: " <<popSize;
	int swapIndex;
	cout << "Rank by index: ";
	for(int i = 0; i < popSize; i++){
		//bestScore = tempIndices.at(i);
		for(int j = i+1; j< popSize; j++){
			if(tempScores.at(i) > tempScores.at(j)){
				swapScore = tempScores.at(i);
				swapIndex = tempIndices.at(i);
				tempScores.at(i) = tempScores.at(j);
				tempIndices.at(i) = tempIndices.at(j);
				tempScores.at(j) = swapScore;
				tempIndices.at(j) = swapIndex; 
			}
		}
		//cout << i << ": ";
		cout << tempIndices.at(i) << " ";
	}
	cout << endl;
	*score = tempScores.at(0);
	*tempScoresSave = Map<const VectorXi, Eigen::Unaligned>(tempScores.data(), tempScores.size());
	//*tempScoresSave = tempScores;
	return tempIndices;
}

void displayPop(vector<crossbar> &population){
	for (int i = 0; i < population.size(); i++){
		cout <<" mask: " << endl << population.at(i).mask() <<endl << endl;
		cout <<" synapses: " << endl << population.at(i).synapses() <<endl;	
	}
}

int activationFunction(int n, int actFunct){ //this is the activation function	
//	if(actFunct == 0){return 0;}
	if(actFunct < 0){
		if(n + actFunct < 0){return 0;}
		else{return n + actFunct;}
	}
	else{
		if(n > actFunct){
			if(n > 15){return 15 - actFunct;}
			return n - actFunct;	
		}
		else{
			return 0;
		}
	}
//	if(n < 7){return n;}
//	else{return 7;}		
}


VectorXi translateInput(crossbar &crossbarArray, VectorXi input){
	
//	while(input.size() < crossbarArray.size()){
//		input.emplace(input.begin(),0);
//	}
	VectorXi output(crossbarArray.size());
	for(int i=0; i < crossbarArray.size(); i++){
		if(i > (crossbarArray.size() - layerBreadth.at(0)) ){ 
			output(i) = input(i - crossbarArray.size() + layerBreadth.at(0));
		}
		else{
			output(i) = 0;
		}
	}
//	cout << output.transpose();
	return output;
}


VectorXi activation(VectorXi resV, VectorXi actVect){ //this takes an input vector and operates an activation function on it
	VectorXi temp(resV.size());
//	cout << "activation on vector size: " << resV.size() <<endl;
	int i = 0;
	while(i < resV.size() - layerBreadth.at(layerBreadth.size() - 1)){
		temp(i) = activationFunction(resV(i), actVect(i));
		i++;
	}
	//while(i < resV.size()){
//		temp(i) = activationFunction(resV(i), false);
//		i++;
//	}
	return temp;	
}

VectorXi forwardProp(crossbar &crossbarArray, VectorXi input){ //this propagates a given input vector forward, returns output
	VectorXi resV = input;
	VectorXi trailingresV = resV;	
	for(int i=0; i < crossbarArray.layerCount() - 1; i++){
	//	cout << resV << endl <<endl;
//		
	//	cout << crossbarArray.synapses() << endl <<endl;		
//cout <<i<< "th iteration started" <<endl;
//		cout <<resV.transpose() << endl;
//		cout << "about to multiply by crossbar" <<endl;
		//cout<< crossbarArray.synapses().cast<int>() * resV;
		trailingresV = crossbarArray.synapses().cast<int>() * resV;
		resV = activation(trailingresV,crossbarArray.actFunct());  
//		cout <<i<<"th iteration done" << endl;
	}
	VectorXi output(layerBreadth[crossbarArray.layerCount() - 1 ]);
	//cout << resV << endl <<endl;

	for(int i=0; i < layerBreadth[crossbarArray.layerCount() -1]; i++){
		output(i) = resV(i);
	} 
//	cout << output.transpose() <<endl;
	return output;
}
	

void initPopulation(vector<crossbar> &population, int size){
	for(int i = 0; i < size; i++){
//		crossbar temp = new crossbar(layerBreadth, 0,1);
	//	population.push_back(temp);
	}	
}

/*
lookupVector impulseLookup(crossbar &crossbarArray, VectorXd resV){
	for(int i=0; i<size; i++){
		resV.setZero();
		resV(i)=1;
		matrixLookup(i) = crossbarArray.synapses() * resV;
			
	}
	
}
*/

MatrixXd reservoirGen(int size, int intercon, int synapseSize){
	cout << "Setting up randomly weighted reservoir: " << size*size << " synapses" << endl;
	srand (time(NULL));
	intercon = 200; //fixed point, is divided by 100 for average 
	MatrixXd reservoir(size,size);
	for (int i=0; i < size; i++){
		for (int j=0; j < size; j++){
			if(rand() % (size*100) <= intercon){
				reservoir(i,j) = rand() % synapse;
			}
		}
	}
	return reservoir;
}

VectorXd randInput(uint8_t period){
	VectorXd temp(size);
	cout << "Setting up random input vector" <<endl;
	for(int i=0; i < size; i++){
		if(rand() % (size) <=size/period){
			temp(i) = 1;
		}
		else{
			temp(i) = 0;
		}
	}
	return temp;
}

/*
void iterate(VectorXd &resV, neuronVector &resN){
	resV(0)=1;
	bool nonzero = true;
	int n = 0;
	int nActive = 0;
	VectorXd trailingResV(size);
	while(nonzero){
		float activation = nActive / size; 
		trailingResV = resV;
		resV.setZero();
		for(int i=0; i<size; i++){
			if(trailingResV(i)){
				resV += ((int) (15*27*100*100)/(size*synapse)) * trailingResV(i)*matrixLookup(i);
			}
			
				
		}
//		resV = ((int) (15*60*100*100)/(size*synapse)) * reservoir * resV;
//		cout << resV.transpose() << endl;
		nonzero = false;
		nActive = 0;	
		for(int i=0; i< size; i++){
			resN(i).add(resV(i));
			//cout << resN.transpose() << std::endl;
			resV(i) = resN(i).buffer();
			nActive += resV(i);
			resN(i).clearBuffer();
			resN(i).update();
		}
		//resV(0)=1;
		if(nActive != 0){nonzero=true;} 
		//system("clear");
		cout << " Reservoir synapses: " << size*size << " " <<n << "th iteration: "<< "| \% active: " << nActive << "/" << size << std::endl; 
		//cout << resN.transpose() << endl;
//		cout << resV.transpose() << endl;
		n++;
	//	delay(1);
	}
	cout << "dead! \n";
}

*/

