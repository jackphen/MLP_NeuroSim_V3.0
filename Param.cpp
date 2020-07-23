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
#include "Param.h"

Param::Param() {

	maxWeight = 1;	// Upper bound of weight value
	minWeight = 0;	// Lower bound of weight value
    /* Optimization method 
    Available option include: "SGD", "Momentum", "Adagrad", "RMSprop" and "Adam" */

	/* SpikUS Mapping parameters */
	redundancyLevel = 1;

	/* InvPowerIt Algorithm parameters */
	numCycles = 6; 
	problemSize = 64;
	arraySize = 64;
	IPIThreshold = 5e-8; 

	/* Hardware parameters */
	arrayWriteType = true; 			  // Ideal or real write for matrix weights (true: real, false: ideal); 	c hardware, false: ideal software)
	numBitInput = 2;       // # of bits of the input data (=1 for black and white data)
	numBitPartialSum = 10;  // # of bits of the digital output (partial weighted sum output)
	numWeightBit = 6;	// # of weight bits (only for pure algorithm, SRAM and digital RRAM hardware)
	numColMuxed = 16;	// How many columns share 1 read circuit (for analog RRAM) or 1 S/A (for digital RRAM)
	numWriteColMuxed = 16;	// How many columns share 1 write column decoder driver (for digital RRAM)
	writeEnergyReport = true;	// Report write energy calculation or not
	NeuroSimDynamicPerformance = true; // Report the dynamic performance (latency and energy) in NeuroSim or not
	relaxArrayCellHeight = 0;	// True: relax the array cell height to standard logic cell height in the synaptic array
	relaxArrayCellWidth = 0;	// True: relax the array cell width to standard logic cell width in the synaptic array
	arrayWireWidth = 32;	// Array wire width (nm)
	processNode = 32;	// Technology node (nm)
	clkFreq = 2e9;		// Clock frequency (Hz)
	numArrays; 			// Number of arrays/synaptic cores. 
}

