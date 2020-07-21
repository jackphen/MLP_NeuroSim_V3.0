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
#include "IO.h"
#include "Train.h"
#include "Test.h"
#include "Mapping.h"
#include "Definition.h"

using namespace std;

int main() {

	printf("Started.\n"); 

	gen.seed(0);

	/* *********************************************************+
	 * 						ARRAY SIZING
	 * **********************************************************/

	int array_type; 
	cout << "Choose array type [1 1T1R, 2 SRAM]: "; cin >> array_type;

	cout << "# of rows: "; 	cin >> arrayIH->arrayRowSize; param->problemSizeRows = arrayIH->arrayRowSize; 
	cout << "# of cols: ";	cin >> arrayIH->arrayColSize; param->problemSizeCols = arrayIH->arrayColSize; 

	if (array_type == 2) {cout << "SRAM precision [bits]: "; cin >> param->numWeightBit; };

	switch(array_type) {
		case 1: arrayIH->Initialization<RealDevice>(); printf("1T1R %dx%d array initialized.\n\n",arrayIH->arrayRowSize,arrayIH->arrayColSize); break; 
		case 2: arrayIH->Initialization<SRAM>(param->numWeightBit); printf("SRAM %dx%d array initialized.\n\n",arrayIH->arrayRowSize,arrayIH->arrayColSize); break; 
	};
	
	/* *********************************************************+
	 * 						LEAKAGE BENCHMARK
	 * **********************************************************/

	/* Initialization of NeuroSim synaptic cores */
	param->relaxArrayCellWidth = 0;
	NeuroSimSubArrayInitialize(subArrayIH, arrayIH, inputParameterIH, techIH, cellIH);

	/* Calculate synaptic core area */
	NeuroSimSubArrayArea(subArrayIH);

	/* Calculate synaptic core standby leakage power */
	NeuroSimSubArrayLeakagePower(subArrayIH);
	
	/* Initialize the neuron peripheries */
	NeuroSimNeuronInitialize(subArrayIH, inputParameterIH, techIH, cellIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
	
	/* Calculate the area and standby leakage power of neuron peripheries below subArrayIH */
	double heightNeuronIH, widthNeuronIH;
	NeuroSimNeuronArea(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH, &heightNeuronIH, &widthNeuronIH);
	double leakageNeuronIH = NeuroSimNeuronLeakagePower(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
	
	/* Print the area of synaptic core and neuron peripheries */
	double totalNeuronAreaIH = adderIH.area + muxIH.area + muxDecoderIH.area + dffIH.area + subtractorIH.area;

	/* *********************************************************+
	 * 						ARRAY PROGRAMMING 
	 * **********************************************************/

	printf("\nStarted array programming phase.\n"); std::string matrix_name; 
	cout << "Insert input file: "; 	cin >> matrix_name;
	cout << "Select ideal or real cell write [1 real, 0 ideal]: "; 	cin >> param->arrayWriteType; 


	/* Initialize weights and map weights to conductances for hardware implementation */
	WeightInitialize(matrix_name);	WeightToConductance();

	/* *********************************************************+
	 * 				EIGENVECTOR ALGORITHM BENCHMARK 
	 * **********************************************************/	 
	Validate();


	/* *********************************************************+
	 * 				       BENCHMARK RESULTS 
	 * **********************************************************/	 
	printf("\nArea, memory array:\t\t\t %.4e mm^2\n", (subArrayIH->usedArea)*1e6);
	printf("Area, peripherals:\t\t\t %.4e mm^2\n", (totalNeuronAreaIH)*1e6);
	printf("Area, total:\t\t\t\t %.4e mm^2\n", (subArrayIH->usedArea + totalNeuronAreaIH)*1e6 );

	printf("\nArea performance:\t\t\t %.4e TOPS/mm^2\n",pow(param->problemSizeCols,2)/((subArrayIH->usedArea + totalNeuronAreaIH)*1e6)/1e12);	

	printf("\nLeakage power, memory array:\t\t %.4e W\n", subArrayIH->leakage);
	printf("Leakage power, peripherals:\t\t %.4e W\n", leakageNeuronIH);
	printf("Leakage power, total:\t\t\t %.4e W\n", subArrayIH->leakage + leakageNeuronIH);
	printf("-----------------------------------------------------------------\n");

	return 0;
}


