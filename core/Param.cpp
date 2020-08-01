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
#include "math.h"
#include "formula.h"
#include "Debug.h"
#include "Param.h"

/* Logger libraries */
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

extern spdlog::logger logger;

Param::Param() {

	maxWeight = 0;	// Upper bound of weight value
	minWeight = INFINITY;	// Lower bound of weight value
    /* Optimization method 
    Available option include: "SGD", "Momentum", "Adagrad", "RMSprop" and "Adam" */

	/* SpikUS Mapping parameters */
	redundancyLevel = 1;

	/* InvPowerIt Algorithm parameters */
	numCycles = 6; 
	problemSize = 64;
	arraySize = 64;
	IPIThreshold = 5e-8; 
	arrayType = 1;
	matrixGains; 
	matrixFiles;
	targetFile; 
	solutionFile; 
	experimentName;
	outputFile; 
	adaptBenchmark = ""; 

	/* Hardware parameters */
	arrayWriteType = true; 			  // Ideal or real write for matrix weights (true: real, false: ideal); 	c hardware, false: ideal software)
	numBitInput = 8;       // # of bits of the input data (=1 for black and white data)
	numBitPartialSum = 10;  // # of bits of the digital output (partial weighted sum output)
	numWeightBit = 6;	// # of weight bits (only for pure algorithm, SRAM and digital RRAM hardware)
	numColMuxed = 1;	// How many columns share 1 read circuit (for analog RRAM) or 1 S/A (for digital RRAM)
	numWriteColMuxed = 1;	// How many columns share 1 write column decoder driver (for digital RRAM)
	writeEnergyReport = true;	// Report write energy calculation or not
	NeuroSimDynamicPerformance = true; // Report the dynamic performance (latency and energy) in NeuroSim or not
	relaxArrayCellHeight = 0;	// True: relax the array cell height to standard logic cell height in the synaptic array
	relaxArrayCellWidth = 0;	// True: relax the array cell width to standard logic cell width in the synaptic array
	arrayWireWidth = 32;	// Array wire width (nm)
	processNode = 32;	// Technology node (nm)
	clkFreq = 2e9;		// Clock frequency (Hz)
	numArrays; 			// Number of arrays/synaptic cores. 
}

void Param::ReadConfigFile(std::string config_file) {

	/* Attempt to read configuration file */
	std::ifstream config (config_file);
	std::string row; 

	if (config.is_open()) {

		logger.info("Configuration file found!");

		while (!config.eof()) {

			getline(config,row); 

			if (row.at(0) == '#') continue; 

			auto rown = explode(row,'=');

			if (rown[0].compare("experimentName") == 0) {
				experimentName = rown[1];
			}
			else if (rown[0].compare("arrayType") == 0) {
				arrayType = stoi(rown[1]);
			}
			else if (rown[0].compare("problemSize") == 0) {
				problemSize = stoi(rown[1]);
			}
			else if (rown[0].compare("numWeightBit") == 0) {
				numWeightBit = stoi(rown[1]);
			}
			else if (rown[0].compare("redundancyLevel") == 0) {	
				redundancyLevel = stoi(rown[1]);
			}
			else if (rown[0].compare("numArrays") == 0) {
				numArrays = stoi(rown[1]);
			}
			else if (rown[0].compare("matrixFiles") == 0) {
				matrixFiles = explode(rown[1],',');
			}
			else if (rown[0].compare("matrixGains") == 0) {
				rown = explode(rown[1],',');
				for(int ii = 0; ii < rown.size(); ii++)  {
					matrixGains.push_back(stod(rown[ii]));
				}
			}
			else if (rown[0].compare("solutionFile") == 0) {
				solutionFile = rown[1];
			}
			else if (rown[0].compare("arrayWriteType") == 0) {
				arrayWriteType = stoi(rown[1]);
			}
			else if (rown[0].compare("numColMuxed") == 0) {
				numColMuxed = stoi(rown[1]);
			}
			else if (rown[0].compare("numCycles") == 0) {
				numCycles = stoi(rown[1]);
			}			
			else if (rown[0].compare("IPIThreshold") == 0) {
				IPIThreshold = stod(rown[1]);
			}
			else if (rown[0].compare("numBitInput") == 0) {
				numBitInput = stoi(rown[1]);
			}
			else if (rown[0].compare("adaptBenchmark") == 0) {
				adaptBenchmark = rown[1];
			}
		}

		LogConfig(); 
	}
	else { 

		logger.warn("No configuration file found!"); 

		/* Experiment name */
		std::cout << "Experiment name:"; std::cin >> experimentName;

		/* Array type */
		std::cout << "Array type [1 1T1R HP, 2 SRAM, 3 1T1R Fairy]: "; std::cin >> arrayType;

		if (arrayType == 2) {std::cout << "SRAM precision [bits]: "; std::cin >> numWeightBit; };


		/* Problem size */
		std::cout << "Problem size: "; 
		std::cin >> problemSize;

		/* Redundancy */
		std::cout << "Redundancy level: "; 
		std::cin >> redundancyLevel; 

		/* Array quantity */
		std::cout << "Number of arrays: "; 
		std::cin >> numArrays;

		/* Matrix program files and gains */
		for (int jj = 0; jj < numArrays; jj++) {
			printf("Array #%d, target matrix: ",jj); 	
			std::cin >> matrixFiles[jj];
			printf("Array #%d, gain: ",jj); 
			std::cin >> matrixGains[jj];
		}

		/* Array write type */
		std::cout << "Cell write type [1 real, 0 ideal]: "; 	
		std::cin >> arrayWriteType;

		/* Save config prompt */
		char save = 'n'; 
		std::cout << "Save configuration for future reuse? [y/n]";
		std::cin >> save;

		if (save == 'y') WriteConfigFile("config.nsim"); 
	}

	/* Output file */
	outputFile = "experiments/" + experimentName + "/breakdown.txt";

	return;
}

void Param::WriteConfigFile(std::string config_file) {
	logger.warn("Missing WriteConfigFile() function. Config file will not be saved.");
}

void Param::LogConfig() {
	logger.info("experimentName: {}",experimentName.c_str());
	logger.info("arrayType: {}",arrayType);
	logger.info("numWeightBit: {}",numWeightBit); 
	logger.info("problemSize: {}",problemSize);
	logger.info("redundancyLevel: {}",redundancyLevel);
	logger.info("numArrays: {}",numArrays);

	for (int jj = 0; jj < matrixFiles.size(); jj++) {
		logger.info("matrixFiles[{}]: {}",jj,matrixFiles[jj]);
	}
	for (int jj = 0; jj < matrixGains.size(); jj++) {
		logger.info("matrixGains[{}]: {}",jj,matrixGains[jj]);
	}
	logger.info("solutionFile: {}",solutionFile);
	logger.info("outputFile: {}",outputFile);
	logger.info("arrayWriteType: {}",arrayWriteType);
	logger.info("numColMuxed: {}",numColMuxed);
	logger.info("IPIThreshold: {}",IPIThreshold);
	logger.info("numCycles: {}",numCycles);
	logger.info("numBitInput: {}",numBitInput);
	logger.info("adaptBenchmark: {}",adaptBenchmark);
}