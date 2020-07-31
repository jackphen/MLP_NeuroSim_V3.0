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

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <vector>
#include "Cell.h"
#include "Array.h"
#include "formula.h"
#include "NeuroSim.h"
#include "Param.h"
#include "Test.h"
#include "Mapping.h"
#include "Debug.h"
#include "Results.h"
#include "Definition.h"

using namespace std;

void InitSystem(); 

char benchmark[200];

int main() {

	TRACE("Started.\n"); 

	/* To be moved in Results::Results() */
	std::vector<double> heightNeurons;		// Vector of heights of neuron peripheries circutry
	std::vector<double> widthNeurons;		// Vector of widths of neuron peripheries circuitry
	std::vector<double> leakageNeurons;		// Vector of leakage powers of neuron peripheries circuitry
	std::vector<double> totalNeuronAreas;	// Vector of total area of neuron peipheries circuitry

	sprintf(benchmark,"benchmark = ["); 

	gen.seed(0);

	/* *********************************************************+
	 * 					SETUP INITIALIZATION
	 * **********************************************************/

	printf("Started\n"); 
	param->ReadConfigFile("config.nsim"); 

	/* Initialize each memory array */
	param->arraySize = (int) sqrt(param->redundancyLevel)*param->problemSize; 

	for (int jj = 0; jj < param->numArrays; jj++) {
		arrays.push_back(new Array(param->arraySize,param->arraySize, param->arrayWireWidth));
		arrays[jj]->arrayRowSize = param->arraySize;
		arrays[jj]->arrayColSize = param->arraySize; 

		switch(param->arrayType) {
			case 1: arrays[jj]->Initialization<PoliMiDevice>(); TRACE("1T1R %dx%d array %d initialized.\n\n",arrays[jj]->arrayRowSize,arrays[jj]->arrayColSize,jj); break; 
			case 2: arrays[jj]->Initialization<SRAM>(param->numWeightBit); TRACE("SRAM %dx%d array %d initialized.\n\n",arrays[jj]->arrayRowSize,arrays[jj]->arrayColSize,jj); break; 
			case 3: arrays[jj]->Initialization<RealDevice>(); TRACE("1T1R %dx%d array %d initialized.\n\n",arrays[jj]->arrayRowSize,arrays[jj]->arrayColSize,jj); break; 
		};
	}


	/* *********************************************************+
	 * 						LEAKAGE BENCHMARK
	 * **********************************************************/

	/* This was fun. If the array is too small, the resulting width
	 * is not enough to accomodate even one subtractor. This leads 
	 * to a divide by zero inside NeuroSim calculateArea() method
	 * for the Subtractor class, which then results in the overflow
	 * of the associated double area attribute. 
	 * 
	 * This bug manifests as a negative benchmarked area for the 
	 * array peripherals.
	 * 
	 * Since we (a) do not want to tinker with NeuroSim excessively,
	 * and (b) consider this to be also a realistic issue, we chose
	 * to relax the array cell width should the array size be too 
	 * small to accomodate one subtractor, and not relax it if at 
	 * least one subtractor can fit along the array. The magic size
	 * of 34 was found by trial and error.
	 * */
	param->relaxArrayCellWidth = (param->arraySize < 34);
	
	/* Perform initialization of each subArray */
	for (int jj = 0; jj < param->numArrays; jj++) {

		subArrays.push_back(new SubArray(inputParameterIH,techIH,cellIH)); 
		adders.push_back(new Adder(inputParameterIH,techIH,cellIH));
		muxs.push_back(new Mux(inputParameterIH,techIH,cellIH));
		muxDecoders.push_back(new RowDecoder(inputParameterIH,techIH,cellIH));
		dffs.push_back(new DFF(inputParameterIH,techIH,cellIH));
		subtractors.push_back(new Subtractor(inputParameterIH,techIH,cellIH));

		/* Initialization of NeuroSim synaptic cores */
		NeuroSimSubArrayInitialize(subArrays[jj], arrays[jj], inputParameterIH, techIH, cellIH);
		TRACE("\tInitialized synaptic core #%d.\n",jj);

		/* Calculate synaptic core area */
		/* This includes the memory array and the "reading" peripherals */
		NeuroSimSubArrayArea(subArrays[jj]);
		TRACE("\tCalculated area of synaptic core #%d.\n",jj);

		/* Calculate synaptic core standby leakage power */
		NeuroSimSubArrayLeakagePower(subArrays[jj]);
		TRACE("\tCalculated leakage of synaptic core #%d.\n",jj);
		
		/* Initialize the neuron peripheries */
		/* These are used for the digital post-processing and feedback operation */
		NeuroSimNeuronInitialize(subArrays[jj], inputParameterIH, techIH, cellIH, *adders[jj], *muxs[jj], *muxDecoders[jj], *dffs[jj], *subtractors[jj]);
		TRACE("\tInitialized peripherals of synaptic core #%d.\n",jj);

		/* Area of neuron peripheries */
		heightNeurons.push_back(0); widthNeurons.push_back(0); 
		NeuroSimNeuronArea(subArrays[jj], *adders[jj], *muxs[jj], *muxDecoders[jj], *dffs[jj], *subtractors[jj], &heightNeurons[jj], &widthNeurons[jj]);
		totalNeuronAreas.push_back((*adders[jj]).area + (*muxs[jj]).area + (*muxDecoders[jj]).area + (*dffs[jj]).area + (*subtractors[jj]).area);
		TRACE("\tCalculated area of perihperals of synaptic core #%d.\n",jj);

		/* Leakage of neuron peripheries */
		leakageNeurons.push_back(NeuroSimNeuronLeakagePower(subArrays[jj], *adders[jj], *muxs[jj], *muxDecoders[jj], *dffs[jj], *subtractors[jj]));
		TRACE("\tCalculated leakage of peripherals of synaptic core #%d.\n",jj);

		TRACE("Synaptic Core #%d initialized.\n\n",jj);
	}

	printf("\nInitialized all arrays.\n"); 
 

	/* *********************************************************+
	 * 						ARRAY PROGRAMMING 
	 * **********************************************************/

	printf("\nStarted array programming phase.\n");

	for(int jj = 0; jj < param->numArrays; jj++) {
		WeightInitialize(param->matrixFiles[jj],arrays[jj]);  
	}

	/* *********************************************************+
	 * 				EIGENVECTOR ALGORITHM BENCHMARK 
	 * **********************************************************/	 
	PowerIt_Full();


	/* *********************************************************+
	 * 				       BENCHMARK RESULTS 
	 * **********************************************************/	 
	printf("-----------------------------------------------------------------\n");
	printf("                           SYSTEM REPORT                         \n");
	printf("-----------------------------------------------------------------\n");

	double totalNeuronArea = 0; 
	double totalMemoryArea = 0; 
	double totalNeuronLeakage = 0; 
	double totalMemoryLeakage = 0; 
	double totalSubArrayArea = 0; 
	for(int jj = 0; jj < param->numArrays; jj++) {
		TRACE("-----------------------------------------------------------------\n");
		TRACE("                           ARRAY #%d                             \n",jj);
		TRACE("-----------------------------------------------------------------\n");
		TRACE("\nArea, synaptic core:\t\t\t %.4e mm^2\n", (subArrays[jj]->usedArea));
		TRACE("\t\tmemory array:\t\t\t %.4e mm^2\n", (subArrays[jj]->areaArray));
		TRACE("\t\tperipherals:\t\t\t %.4e mm^2\n", (subArrays[jj]->usedArea-subArrays[jj]->areaArray));
		TRACE("Area, neuron:\t\t\t %.4e mm^2\n", (totalNeuronAreas[jj])*1e6);
		TRACE("Area, total:\t\t\t\t %.4e mm^2\n", (subArrays[jj]->areaArray + totalNeuronAreas[jj])*1e6 );

		TRACE("\nLeakage power, synaptic core:\t\t %.4e W\n", subArrays[jj]->leakage);
		TRACE("Leakage power, neuron:\t\t\t %.4e W\n", leakageNeurons[jj]);
		TRACE("Leakage power, total:\t\t\t %.4e W\n", subArrays[jj]->leakage + leakageNeurons[jj]);
		TRACE("-----------------------------------------------------------------\n");

		totalNeuronArea += totalNeuronAreas[jj];
		totalMemoryArea += subArrays[jj]->areaArray;
		totalNeuronLeakage += leakageNeurons[jj];
		totalMemoryLeakage += subArrays[jj]->leakage; 
		totalSubArrayArea += subArrays[jj]->usedArea; 
	}

	TRACE("-----------------------------------------------------------------\n");
	TRACE("                               TOTAL                             \n");
	TRACE("-----------------------------------------------------------------\n");
	printf("\nArea, synaptic cores:\t\t\t %.4e mm^2\n", (totalSubArrayArea)*1e6);
	printf("\tmemory arrays:\t\t\t\t %.4e mm^2\n", (totalMemoryArea)*1e6);
	printf("\tperipherals:\t\t\t\t %.4e mm^2\n", (totalSubArrayArea - totalMemoryArea)*1e6);
	printf("Area, neuron:\t\t\t\t %.4e mm^2\n", (totalNeuronArea)*1e6);
	printf("Area, total:\t\t\t\t %.4e mm^2\n", (totalSubArrayArea + totalNeuronArea)*1e6 );

	printf("\nPerformance Density:\t\t\t %.4e TOPS/mm^2\n",results->Throughput/((totalSubArrayArea + totalNeuronArea)*1e6)/1e12);	

	printf("\nLeakage power, synaptic cores:\t\t %.4e W\n", totalMemoryLeakage);
	printf("Leakage power, neurons:\t\t\t %.4e W\n", totalNeuronLeakage);
	printf("Leakage power, total:\t\t\t %.4e W\n", totalMemoryLeakage + totalNeuronLeakage);
	printf("-----------------------------------------------------------------\n");

	sprintf(benchmark,"%s%.4e,",benchmark,totalMemoryArea); 
	sprintf(benchmark,"%s%.4e,",benchmark,totalNeuronArea); 
	sprintf(benchmark,"%s%.4e,",benchmark,totalMemoryArea+totalNeuronArea); 
	sprintf(benchmark,"%s%.4e,",benchmark,totalMemoryLeakage); 
	sprintf(benchmark,"%s%.4e,",benchmark,totalNeuronLeakage); 
	sprintf(benchmark,"%s%.4e];\n",benchmark,totalMemoryLeakage+totalNeuronLeakage); 

	printf("\n\n%s",benchmark);
	TRACE("\nCompleted!\n");

	results->PerformanceDensity = results->Throughput/(totalSubArrayArea + totalNeuronArea);
	results->TotalArea = totalSubArrayArea + totalNeuronArea;
	results->TotalCoreArea = totalSubArrayArea; 
	results->TotalCoreMemoryArea = totalMemoryArea; 
	results->TotalCorePeripheralArea = totalSubArrayArea - totalMemoryArea; 



	for (int jj = 0; jj < param->numArrays; jj++) {
		(*subArrays[jj]).PrintProperty(); 
	}

	return 0;
}