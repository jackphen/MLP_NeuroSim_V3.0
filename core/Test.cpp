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
#include "formula.h"
#include "Param.h"
#include "Array.h"
#include "Mapping.h"
#include "NeuroSim.h"
#include "Cell.h"
#include "Debug.h"
#include "Results.h"

/* Logger libraries */
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

extern spdlog::logger logger;

extern Param *param;
extern Results *results; 

extern Technology techIH;
extern Array *arrayIH;
extern SubArray *subArrayIH;
extern Adder adderIH;
extern Mux muxIH;
extern RowDecoder muxDecoderIH;
extern DFF dffIH;
extern Subtractor subtractorIH;

extern std::vector<Array* > arrays;
extern std::vector<SubArray* > subArrays;
extern std::vector<Adder*> adders;
extern std::vector<Mux*> muxs;
extern std::vector<RowDecoder*> muxDecoders;
extern std::vector<DFF*> dffs;
extern std::vector<Subtractor*> subtractors;

extern char benchmark[200]; 

/* Validation */
void PowerIt_Full() {

	int numBatchReadSynapse;    // # of read synapses in a batch read operation (decide later)
	
	std::vector<std::vector<double>> HWIout(param->numArrays,
										  std::vector<double>(arrays[0]->arrayColSize,0));		// Column output current vector
	std::vector<std::vector<int>> 	 HWInput(param->numArrays,
										  std::vector<int>(arrays[0]->arrayRowSize,0));		// Hardware input vector (after redunancy reshape)
	std::vector<std::vector<double>> HWOutput(param->numArrays,
										  std::vector<double>(arrays[0]->arrayColSize,0));		// Hardware output vector (before redunancy reshape)

	std::vector<double> temp_Iout(arrays[0]->arrayColSize,0); 		// Temporary vector to store slicing rebuild

	std::vector<double> Input(param->problemSize,0); 
	std::vector<double> Output(param->problemSize,0); 
	std::vector<double> EvExact(param->problemSize,0); 					
	
	std::vector<double> sumArrayReadEnergy(param->numArrays,0); 		// Temporary variable vector to store array read energy for each subarray
	std::vector<double> sumSubArrayReadEnergy(param->numArrays,0); 		// Temporary variable vector to store subarray (memory array + readout peripherals) read energy for each subarray
	std::vector<double> sumNeuroSimReadEnergy(param->numArrays,0);		// Temporary variable vector to store neuron peripherals read energy for each subarray
	std::vector<double> sumReadLatency(param->numArrays,0); 			// Temporary variable vector to store overall latency for each subarray

	double maxIout = 0; 		// Maximum output current in a single iteration (needed to quantize)
	double maxOutput = 1; 		// Maximum output value (needed to normalize)
	double max_ev = 0; 			// Maximum value of the eigenvector (needed to normalize)

    double readVoltageIH;		// Cell read voltage (needed for energy computation)				
    double readPulseWidthIH;	// Cell read pulse width (needed for energy computation)
    
	
	/* All arrays have the same memory device, so readVoltage and readPulse are read from the first. */
	if(eNVM* temp = dynamic_cast<eNVM*>(arrays[0]->cell[0][0]))
    {
        readVoltageIH = static_cast<eNVM*>(arrays[0]->cell[0][0])->readVoltage;
        readPulseWidthIH = static_cast<eNVM*>(arrays[0]->cell[0][0])->readPulseWidth;
    }

	logger.info("Initiating eigenvector algorithm.");

	/* Threshold adjustment: the nominal threshold is adjusted to the mapping conductance */
	bool digitalNVM = false, parallelRead = false;

	if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrays[0]->cell[0][0])) {
		logger.trace("Threshold nominal: {}, Gmax: {}", param->IPIThreshold,(static_cast<eNVM*>(arrays[0]->cell[0][0])->maxConductance)*readVoltageIH);
		param->IPIThreshold = param->IPIThreshold * readVoltageIH / param->maxWeight * (static_cast<eNVM*>(arrays[0]->cell[0][0])->maxConductance);
		logger.trace("Threshold real: {}",param->IPIThreshold);
	}
	else if(DigitalNVM*temp = dynamic_cast<DigitalNVM*>(arrays[0]->cell[0][0])) {
		
		digitalNVM = true;
		if(static_cast<DigitalNVM*>(arrays[0]->cell[0][0])->parallelRead == true) parallelRead = true;

		if(digitalNVM && parallelRead) {} // TODO: digitalNVM
		else { // SRAM
			param->IPIThreshold = param->IPIThreshold / param->maxWeight * (pow(2,param->numWeightBit)-1);
		}
	}		

	/* Read exact solution from file */
	std::ifstream eigenvector("experiments/" + param->experimentName + "/" + param->solutionFile);
    std::string row; 

	getline(eigenvector,row); 
	auto rown = explode(row,';');

	for (int j = 0; j < param->problemSize; j++) {
		EvExact[j] = stof(rown[j]); 
		max_ev = MAX(abs(EvExact[j]),max_ev); 
	}
	for (int j = 0; j < param->problemSize; j++) {
		EvExact[j] = EvExact[j] / max_ev; 
	}

    eigenvector.close();

	for (int ii = 0; ii < param->numArrays; ii++) {
		/* ****************************************************************
		* Step 1: row/gate precharge energy. 
		* This step is needed to compute the energy used to charge 
		* either the gate capacitances (1T1R) or the TE lines (X-point)
		* of unselected cells. 
		* ****************************************************************/
		if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrays[ii]->cell[0][0])) {  // Analog eNVM
			if (static_cast<eNVM*>(arrays[ii]->cell[0][0])->cmosAccess) {  // 1T1R
				sumArrayReadEnergy[ii] += arrays[ii]->wireGateCapRow * techIH.vdd * techIH.vdd * (arrays[ii]->arrayColSize*arrays[ii]->arrayRowSize); // All WLs open
			}
		} else if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrays[ii]->cell[0][0])) { // Digital eNVM
			if (static_cast<eNVM*>(arrays[ii]->cell[0][0])->cmosAccess) {  // 1T1R
				sumArrayReadEnergy[ii] += arrays[ii]->wireGateCapRow * techIH.vdd * techIH.vdd;  // Selected WL
			} else {    // Cross-point
				sumArrayReadEnergy[ii] += arrays[ii]->wireCapRow * techIH.vdd * techIH.vdd * (arrays[ii]->arrayColSize - 1);    // Unselected WLs
			}
		}
	}

	/* Inverse Power Iteration*/
	for (int i = 0; i < param->numCycles; i++)
	{
		/* Integerize and create input vector for each array. */
		for (int ii = 0; ii < param->numArrays; ii++)
		{
			for (int r = 0; r < (int) sqrt(param->redundancyLevel); r++) 
			{
				for (int j = 0; j < param->problemSize; j++) 
				{
					HWInput[ii][j+r*param->problemSize] = integerize(Input[j],param->numBitInput);
					logger.trace("Input[{}] = {} <--> HWInput[{}+{}] = {}",j,Input[j],j,r*param->problemSize,HWInput[ii][j+r*param->problemSize]);
				}
			}
		}

		/* Perform Parallel Matrix-Vector Multiplication (PMVM). */
		for (int ii = 0; ii < param->numArrays; ii++) 
		{

			/* Each MVM is computed column by column, row by row. */
			for (int j = 0; j < arrays[ii]->arrayColSize; j++) {
			
				/* ****************************************************************
				* Step 2: MVM. 
				* The input vector is fed to the network bit-by-bit and the output
				* current is computed along with the read energy on the selected
				* column. 
				* ****************************************************************/
				for (int n = 0; n < param->numBitInput; n++) {

					/* AnalogNVM */
					if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrays[ii]->cell[0][0])) { 
						
						double Isum = 0;    // weighted sum current
						double Isum_logic = 0; 		// Weighted sum current, considering shift&add					

						// For each row, calculate current contribution accumulating on the column. 
						for (int k = 0; k < arrays[ii]->arrayRowSize; k++) {
							if ((HWInput[ii][k] >> n) & 1) {    // If the nth bit of the input on row k is 1, accumulate current
								Isum += (arrays[ii]->ReadCell(j,k));
								Isum_logic += (arrays[ii]->ReadCell(j,k))*pow(2,n); 
								sumArrayReadEnergy[ii] += arrays[ii]->wireCapRow * readVoltageIH * readVoltageIH;   // Selected BLs (1T1R) or Selected WLs (cross-point)
							}
						}

						// Reading this bit required energy. 
						sumArrayReadEnergy[ii] += Isum * readVoltageIH * readPulseWidthIH;

						// Once energy calculation is complete, the value of Isum we care about for algorithm purposes is stored in Isum_logic. 
						Isum = Isum_logic; 

						// The computed current is the one corresponding to column j, hence to output j.
						HWIout[ii][j] = Isum;
					}
					else { 
						if(digitalNVM && parallelRead) // parallel read-out for DigitalNVM
						{
							/* TODO: insert algorithm using digitalNVM and parallel readout */
						}
						else {	 // Digital NVM or SRAM row-by-row readout				
							int Dsum = 0;
							int Dsum_logic = 0; 

							for (int k=0; k < arrays[ii]->arrayRowSize; k++) {
								if ((HWInput[ii][k] >> n) & 1) {    // if the nth bit of HWInput[k] is 1
									Dsum += (int)(arrays[ii]->ReadCell(j,k));
									Dsum_logic += (int)(arrays[ii]->ReadCell(j,k))*pow(2,n);
								}
							}
							if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrays[ii]->cell[0][0])) {    // Digital eNVM
								sumArrayReadEnergy[ii]  += static_cast<DigitalNVM*>(arrays[ii]->cell[0][0])->readEnergy * arrays[ii]->numCellPerSynapse * arrays[ii]->arrayRowSize;
							} 
							else {    // SRAM
								sumArrayReadEnergy[ii]  += static_cast<SRAM*>(arrays[ii]->cell[0][0])->readEnergy * arrays[ii]->numCellPerSynapse * arrays[ii]->arrayRowSize;
							}

							// Once energy calculation is complete, the value of Dsum we care about for algorithm purposes is stored in Dsum_logic. 
							Dsum = Dsum_logic; 

							// The computed current is the one corresponding to column j, hence to output j.
							HWIout[ii][j] = (double) (Dsum);//*2.1790e-4/4944;
						}
					}
				}
			}
		}

		logger.trace("MVM Complete. Initiating output extraction."); 

		/* ********************************************
		 * 				Extraction of Output
		 * 
		 * The output is quantized with numBitPartialSum
		 * precision, and using the maximum current read
		 * as full-scale reference. 
		 * ********************************************/

		/* First, currents are properly combined through gains. */
		
		for (int ii = 0; ii < param->numArrays; ii++) {
			for (int j = 0; j < arrays[ii]->arrayColSize; j++) {
				temp_Iout[j] += HWIout[ii][j]*param->matrixGains[ii];
			}
		}

		/* Then, redundancy is folded and averaged. */
		for (int r = 0; r < (int) sqrt(param->redundancyLevel); r++) {
			for (int j = 0; j < param->problemSize; j++) {
				Output[j] += temp_Iout[j+r*param->problemSize];
			}
		}
		for (int j = 0; j < param->problemSize; j++) {
			Output[j] = Output[j] / param->redundancyLevel;
			logger.trace("Output[{}] = {}",j,Output[j]); 
		}

		/* Output is then scaled by Ith. */
		maxOutput = 0;
		for (int j = 0; j < param->problemSize; j++) {
			Output[j] = Output[j] / param->IPIThreshold; 
			if (Output[j] > maxOutput) maxOutput = Output[j];
		}

		/* Output is normalized wrt maximum value. */
		if (maxOutput == 0) maxOutput = 1;
		if (maxOutput > 1) {
			for (int j = 0; j < param->problemSize; j++) {
				Output[j] = Output[j] / maxOutput; 
			}
		}

		/* Normalized output is then quantized (sampling by ADC). */
		for (int j = 0; j < param->problemSize; j++) {
			Output[j] = quantize(Output[j],param->numBitPartialSum); 
		}

		
		/* ********************************************
		 * 					Feedback
		 * 
		 * Output in this time step is re-applied to
		 * the network in the next time step. 
		 * 
		 * If we are at the initial time step, a delta 
		 * input is provided. 
		 * ********************************************/

		for (int j = 0; j < param->problemSize; j++) {
			Input[j] = Output[j] + (i == 0);
		}

		/* ********************************************
		 * 					Report
		 * 
		 * After each iteration, we compute the accuracy
		 * of the solution and update the energy report.
		 * ********************************************/

		// Mean absolute error: results->MAE = sum(|x-x_hat|)/N
		results->MAE = 0; 
		results->MRE = 0; 
		for (int j = 0; j < param->problemSize; j++) {
			results->MAE += abs(Output[j] - EvExact[j]);
			results->MRE += abs(Output[j] - EvExact[j]) / abs(EvExact[j]);
		}
		results->MAE = results->MAE / param->problemSize;
		results->MRE = results->MRE / param->problemSize; 

		// Energy report
		for (int ii = 0; ii < param->numArrays; ii++) {
			numBatchReadSynapse = (int)ceil((double) arrays[ii]->arrayColSize/param->numColMuxed);
			#pragma omp critical    // Use critical here since NeuroSim class functions may update its member variables
			for (int j=0; j<arrays[ii]->arrayColSize; j+=numBatchReadSynapse) {
				int numActiveRows = 0;  // Number of selected rows for NeuroSim
				for (int n=0; n<param->numBitInput; n++) {
					for (int k=0; k<arrays[ii]->arrayRowSize; k++) {
						if ((HWInput[ii][k] >> n) & 1) {    // if the nth bit of HWInput on row k is 1
							numActiveRows++;
						}
					}
				}

				subArrays[ii]->activityRowRead = (double)numActiveRows/arrays[ii]->arrayRowSize/param->numBitInput;
				sumSubArrayReadEnergy[ii] += NeuroSimSubArrayReadEnergy(subArrays[ii]);
				sumNeuroSimReadEnergy[ii] += NeuroSimNeuronReadEnergy(subArrays[ii], *adders[ii], *muxs[ii], *muxDecoders[ii], *dffs[ii], *subtractors[ii]);
				sumReadLatency[ii] += NeuroSimSubArrayReadLatency(subArrays[ii]);
				sumReadLatency[ii] += NeuroSimNeuronReadLatency(subArrays[ii], *adders[ii], *muxs[ii], *muxDecoders[ii], *dffs[ii], *subtractors[ii]);
			}

			logger.trace("Read latency[{}]: {} s", ii, sumReadLatency[ii]);
			logger.trace("Read energy[{}]: {} J",ii, sumNeuroSimReadEnergy[ii]);
		}

		logger.trace("Accuracy at cycle #{}: {}", i, results->MAE);
	}	


	logger.info("Eigenvector algorithm completed successfully.");
	
	std::string vcomp = "v_comp = ["; 
	for(int j = 0; j < param->problemSize; j++) {
		vcomp = vcomp + std::to_string(Output[j]) + ",";
	}
	vcomp = vcomp + "]';";
	logger.warn("{}",vcomp);

	/* *********************************************************+
	 * 				       BENCHMARK RESULTS 
	 * **********************************************************/	 
	
	for(int ii = 0; ii < param->numArrays; ii++) {
		if (sumReadLatency[ii] > results->totalLatency) results->totalLatency = sumReadLatency[ii];
		results->totalCoreMemoryEnergy += sumArrayReadEnergy[ii];
		results->totalCorePeripheralEnergy += sumSubArrayReadEnergy[ii];
		results->totalCoreEnergy += (sumArrayReadEnergy[ii] + sumSubArrayReadEnergy[ii]);
		results->totalNeuronEnergy += sumNeuroSimReadEnergy[ii];	
	}
	results->totalEnergy = results->totalCoreEnergy + results->totalNeuronEnergy; 
	
	return; 
}