/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Coumons Attribution-NonCoumercial 4.0 International Public License 
* http://creativecoumons.org/licenses/by-nc/4.0/legalcode.
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

#ifndef RESULTS_H_
#define RESULTS_H_

/** Class to hold benchmark results and generate export files. */
class Results {
public:
	Results();
	void Benchmark(); 
	void GenerateCSV();
	void PrintResults(); 
	void PrintResultsToFile(const char* outstr);

	double MAE; 				// Mean Absolute Error of the computed solution
	double totalFlops; 			// Total number of floating-point operations

	double throughput;			// Measured throughput in FLOPS 
	double energyPerformance; 	// Measured energy performance in FLOPS/W
	double performanceDensity;  // Measured performance density in FLOPS/um^2
	double SWaP; 				// Measured SWaP (Space-Wattage Performance) in FLOPS/(W um^2)
	
	double totalEnergy; 		// Total measured energy of the system in J
	double totalCoreEnergy; 	// Total measured energy of the synaptic cores in J
	double totalNeuronEnergy;	// Total measured energy of the neurons in J
	double totalCorePeripheralEnergy; 	// Total measured energy of synaptic cores' peripherals in J
	double totalCoreMemoryEnergy; 		// Total measured energy of synaptic cores' memory arrays in J

	double totalLatency; 		// Total latency in s

	double totalArea; 			// Total occupied area in um^2
	double totalCoreArea; 		// Total area of the synaptic cores in um^2
	double totalNeuronArea; 	// Total area of the neurons in um^2
	double totalCorePeripheralArea; 	// Total area of the synaptic cores' peripherals in um^2
	double totalCoreMemoryArea; 		// Total area of the synaptic cores' memory arrays in um^2

	double totalNeuronLeakage; 
	double totalCoreLeakage; 
};

#endif
