/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
*   
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer. 
*   
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
*   
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <string>
#include <sstream>
#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <vector>


#ifndef PARAM_H_
#define PARAM_H_

class Param {
public:
	Param();
	void ReadConfigFile(std::string config_file); 
	void WriteConfigFile(std::string config_file); 
	void LogConfig();

	/* SpikUS Mapping parameters */
	int redundancyLevel; 	// 1: single copy, 4: 4x copy of a single matrix in one array ([A1 A2; A3 A4])
	int arraySize;
	double maxWeight;	// Upper bound of weight value
	double minWeight;	// Lower bound of weight value
    

	/* InvPowerIt Algorithm parameters */
	int numCycles;  		// # of network cycles
	int problemSize; 
	double IPIThreshold; 
	int arrayType; 
	std::vector<double> matrixGains;

	std::vector<std::string> matrixFiles;
	std::string targetFile; 
	std::string solutionFile; 
	std::string experimentName; 
	std::string outputFile; 

	double throughput; 
	double pdensity; 
	double eperformance; 

	/* Hardware parameters */
	bool arrayWriteType; 			// Ideal or real write for matrix weights (true: real, false: ideal); 
	int numBitInput;		// # of bits of the input data (=1 for black and white data)
	int numBitPartialSum;	// # of bits of the digital output (partial weighted sum output)
	int numWeightBit;	// # of weight bits (only for pure algorithm, SRAM and digital RRAM hardware)
	int numColMuxed;	// How many columns share 1 read circuit (for analog RRAM) or 1 S/A (for digital RRAM)
	int numWriteColMuxed;	// How many columns share 1 write column decoder driver (for digital RRAM)
	bool writeEnergyReport;	// Report write energy calculation or not
	bool NeuroSimDynamicPerformance; // Report the dynamic performance (latency and energy) in NeuroSim or not
	bool relaxArrayCellHeight;	// True: relax the array cell height to standard logic cell height in the synaptic array
	bool relaxArrayCellWidth;	// True: relax the array cell width to standard logic cell width in the synaptic array
	double arrayWireWidth;	// Array wire width (nm)
	int processNode;	// Technology node (nm)
	double clkFreq;		// Clock frequency (Hz)
	int numArrays; 		// Number of separate arrays/synaptic cores
};

#endif
