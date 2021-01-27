// This is a small demonstration program showing how the Card and Deck classes are used.
#include <iostream>
#include <cstdint>
//#include <cstdlib>
#include <string>
#include <Eigen/Dense>
//#include "spikeneuron.h" 
#include <time.h>
#include <vector>
#include "crossbar.h"
#include <stdlib.h>
#include "trainingdata.h"

using namespace std;
using namespace Eigen;


// PROTOTYPES for functions used by this demonstration program:
TrainingData::TrainingData(string fileName, int startcol, int inputsize, int outputsize){
	ifstream myFile(fileName);
	string line;
	vector<string> chunk;	
	string word;
	int lineIndex = 0;
	char currentChar;
	getline(myFile, line);
	int delimiterIndex = 0;
	bool dataFlag = false;
	int outputRow;
	int inputRow;
	int charIndex;
	uint8_t tempChar;
	int d;
	while (getline(myFile, line)) {
//		cout << line << endl;
		VectorXi output(outputsize);
		VectorXi input(inputsize);
		outputRow = 0;
		inputRow = 0;
		charIndex = 0;
		while(charIndex < line.size()){
			d = delimiterIndex - startcol + 1;
			tempChar = line[charIndex];
			if(tempChar != ',' && d >= 0 && d < inputsize + outputsize){
			//if(tempChar != ',' && inputRow + outputRow < inputsize + outputsize){
				cout << "num = " << tempChar << " d = " << d << " inputrow, outputrow "  << inputRow << " " << outputRow << endl;
				if(d < inputsize){
					//cout << inputRow << endl;
					input(inputRow) = tempChar - '0';
					inputRow++;
				}	
				else{
					//cout << outputRow << endl;
					output(outputRow) = tempChar - '0';
					outputRow++;					
				}
			}	
			else if(tempChar == ','){
				delimiterIndex++;
			}
			charIndex++;
			/*
			if(delimiterIndex >= startcol && dataFlag){
				if(delimiterIndex <= startcol + inputsize + outputsize){
//					cout << tempChar << " [" << delimiterIndex << "," << startcol + inputsize + outputsize << "]    ";
					if(delimiterIndex > startcol + inputsize){
						output(outputRow) = tempChar - '0';
						//output(outputRow) = 32;
//						cout << output(outputRow) << "   " << outputRow <<endl;
						//cout << tempChar << " ";
						outputRow++;
					}
					else{
						input(inputRow) = tempChar - '0';
						//cout << tempChar << " ";
						inputRow++;
					}
				}
				dataFlag = false;
			} 
			if(tempChar == ','){
				delimiterIndex++; 
				dataFlag = true;
				//cout << "delim ";
				//if(delimiterIndex < startcol + inputsize + outputsize){break;}
			}
			charIndex++;
			*/
		}
		cout << input.transpose() << "->" << output.transpose(); 
		cout << endl;
		delimiterIndex = 0;
		//cout << input.transpose() << "  " << output.transpose() << endl;
		inputs.push_back(input);
		outputs.push_back(output);
		lineIndex++;
	}
	lineIndex = 0;
	numSamples = inputs.size();
	myFile.close();	
}

