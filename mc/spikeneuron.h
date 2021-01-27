//File Name: card.h
//
//Written by Owen Astrachan and Roger Priebe
// This class represents a playing card, i.e., "ace of spades"
// a Card is constructed from a rank (int in range 1..13)
// and a suit (Card::spades, Card::hearts, Card::diamonds,
// Card::clubs)
//
// Cards should be created by a Deck (see deck.h), a Deck returns
// good cards
// The function toString() converts a card to a string, e.g., to print
//
// Accessor functions include
//
// int GetRank()      -- returns 1, 2, ..., 13 for ace, two, ..., king
//
// bool SameSuitAs(c) -- returns true if same suit as Card c
//
// string suitString() -- returns "s", "h", "d" or "c"
//
// Note that the Ace is represented by 1 and the King by 13

#ifndef _SNEURON_H
#define _SNEURON_H

#include <iostream>
#include <cstdint>
#include <string>
using namespace std;

class spikeNeuron
{
  public:
    
	spikeNeuron(){ //default constructor
		threshold = 1000;
		refractDecay = -200;
		leakage = -10;
		outputBuffer = 0;
		integrated = 0;
		refract = false;
		
	}
	spikeNeuron(int _threshold, int _refractDecay, int _leakage){ //custom constructor
		threshold = _threshold;
		refractDecay = _refractDecay;
		leakage = _leakage;
		outputBuffer = 0;
		integrated = 0;
		refract = false;
	}
	void add(int added){ //adds said value to be integrated by the neuron
		if(refract==false){
			integrated += added;
			if(integrated > threshold){
				outputBuffer +=1;
				refract = true;
//				cout << "f";
			}	
		}	
	}	
	int buffer(){ //gives the number of firings before clear
		return outputBuffer;
	}
	void clearBuffer(){ //clears number of times fired
		outputBuffer = 0;
	}
	void update(){ //increments a timestep forward for the neuron
		if(refract){
			if(integrated<1){ //if the neuron has reset, then put it back into integrate mode
				integrated = 0;
				refract = false;
			}
		//	outputBuffer += 1;
			integrated += refractDecay;
		}else{
			integrated += leakage;	
			if(integrated < 1){integrated = 0;}
		}
	}
	//ostream& operator<<(ostream& out, const spikeNeuron& s);
//	friend ostream& operator<<(ostream& out, const spikeNeuron& s){ //display
//		cout << s.outputBuffer;
//		//cout << "I: " << integrated << "/" << threshold << " ";
//	}
	uint8_t operator +(uint8_t offset){ 
		add(offset);
		return buffer();		
	}
    
  private:
	int threshold;
	int refractDecay;
	int integrated;	
	int leakage;
	bool refract;
    	int outputBuffer;
};




#endif
