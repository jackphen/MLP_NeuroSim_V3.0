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
#include <vector>
#include <fstream>
#include <random>
#include "formula.h"
#include "Param.h"
#include "Array.h"
#include "Mapping.h"
#include "NeuroSim.h"
#include "Cell.h"

extern Param *param;

extern std::vector< std::vector<double> > testInput;
extern std::vector< std::vector<int> > dTestInput;
extern std::vector< std::vector<double> > testOutput;

extern std::vector< std::vector<double> > weight1;
extern std::vector< std::vector<double> > weight2;

extern Technology techIH;
extern Array *arrayIH;
extern SubArray *subArrayIH;
extern Adder adderIH;
extern Mux muxIH;
extern RowDecoder muxDecoderIH;
extern DFF dffIH;
extern Subtractor subtractorIH;

extern int correct;		// # of correct prediction

/* Validation */
void Validate() {
	int numBatchReadSynapse;    // # of read synapses in a batch read operation (decide later)
	
	double Iout[param->problemSizeCols]; 		// Column output current vector
	double Input[param->problemSizeCols];		// Input digit vector
	int    HWInput[param->problemSizeCols];		// Hardware input vector (considering redundancy/error arrays)
	double Output[param->problemSizeCols];		// Output digit vector
	double maxIout = 0; 						// Maximum output current in a single iteration (needed to quantize)
	double Ith = param->IPIThreshold; 			// Inverse power iteration threshold
	double ev_exact[param->problemSizeCols]; 	// Exact eigenvector
	double MAE = 0; 							// Mean absolute error

	double sumArrayReadEnergyIH = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
	double sumNeuroSimReadEnergyIH = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
	double sumReadLatencyIH = 0;    // Use a temporary variable here since OpenMP does not support reduction on class member
    double readVoltageIH,readVoltageMSB;
    double readPulseWidthIH,readPulseWidthMSB;
    
	if(eNVM* temp = dynamic_cast<eNVM*>(arrayIH->cell[0][0]))
    {
        readVoltageIH = static_cast<eNVM*>(arrayIH->cell[0][0])->readVoltage;
        readPulseWidthIH = static_cast<eNVM*>(arrayIH->cell[0][0])->readPulseWidth;
		readVoltageIH = 0.05; 
    }
    
	printf("Initiating eigenvector algorithm.\n");
	for (int j = 0; j < param->problemSizeCols; j++) {
		Input[j] = 0; 
		Output[j] = 0; 
	}

	// Read exact solution from file
	std::ifstream eigenvector ("e_exact.txt");
    std::string row; 
	getline(eigenvector,row); 
	auto rown = explode(row,';');

	double max_ev = 0; 
	for (int j = 0; j < param->problemSizeCols; j++) {
		ev_exact[j] = stof(rown[j]); 
		if (abs(ev_exact[j] > max_ev)) max_ev = abs(ev_exact[j]);
	}
	for (int j = 0; j < param->problemSizeCols; j++) {
		ev_exact[j] = ev_exact[j] / max_ev; 
	}

    eigenvector.close();

	for (int i = 0; i < param->numCycles; i++)
	{

		printf("\t Iteration #%d:\n",i); 
		
		/* *******************************************************
		 * 				Matrix Vector Multiplication
		 * 
		 * 
		 * MVM is computed column by column. 
		 * Hence, we sweep along the columns 
		 * and compute the corresponding 
		 * output current and readout energy. 
		 * 
		 * Since we are in the same time-step, 
		 * the same input is applied every time
		 * we compute the current on each column.
		 * 
		 * Hence, the input is binarized at the
		 * beginning of each algorithm iteration
		 * to preserve the speed of the program. 
		 * *******************************************************/

		for (int j = 0; j < param->problemSizeCols; j++) {
			HWInput[j] = integerize(Input[j],param->numBitInput); 
		}

		for (int j=0; j<param->problemSizeCols; j++) {
			
			/* ****************************************************************
			 * Step 1: row/gate precharge energy. 
			 * This step is needed to compute the energy used to charge 
			 * either the gate capacitances (1T1R) or the TE lines (X-point)
			 * of unselected cells. 
			 * ****************************************************************/
			if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayIH->cell[0][0])) {  // Analog eNVM
				if (static_cast<eNVM*>(arrayIH->cell[0][0])->cmosAccess) {  // 1T1R
					sumArrayReadEnergyIH += arrayIH->wireGateCapRow * techIH.vdd * techIH.vdd * param->nInput; // All WLs open
				}
			} else if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrayIH->cell[0][0])) { // Digital eNVM
				if (static_cast<eNVM*>(arrayIH->cell[0][0])->cmosAccess) {  // 1T1R
					sumArrayReadEnergyIH += arrayIH->wireGateCapRow * techIH.vdd * techIH.vdd;  // Selected WL
				} else {    // Cross-point
					sumArrayReadEnergyIH += arrayIH->wireCapRow * techIH.vdd * techIH.vdd * (param->nInput - 1);    // Unselected WLs
				}
			}
		
			/* ****************************************************************
			 * Step 2: MVM. 
			 * The input vector is fed to the network bit-by-bit and the output
			 * current is computed along with the read energy on the selected
			 * column. 
			 * ****************************************************************/
			for (int n = 0; n < param->numBitInput; n++) {

				/* AnalogNVM */
				if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(arrayIH->cell[0][0])) { 
					
					double Isum = 0;    // weighted sum current
					double Isum_logic = 0; 		// Weighted sum current, considering shift&add					

					// For each row, calculate current contribution accumulating on the column. 
					for (int k=0; k<param->problemSizeRows; k++) {

						if ((HWInput[k]>>n) & 1) {    // If the nth bit of the input on row k is 1, accumulate current
							Isum += (arrayIH->ReadCell(j,k));
							Isum_logic += (arrayIH->ReadCell(j,k))*pow(2,n); 
							sumArrayReadEnergyIH += arrayIH->wireCapRow * readVoltageIH * readVoltageIH;   // Selected BLs (1T1R) or Selected WLs (cross-point)
						}
					}

					// Reading this bit required energy. 
					sumArrayReadEnergyIH += Isum * readVoltageIH * readPulseWidthIH;

					// Once energy calculation is complete, the value of Isum we care about for algorithm purposes is stored in Isum_logic. 
					Isum = Isum_logic; 

					// The computed current is the one corresponding to column j, hence to output j.
					Iout[j] = (Isum);

					// Update maximum column current if needed. 
					if (Iout[j] > maxIout) maxIout = Iout[j];
				}
				else { 
					bool digitalNVM = false; 
					bool parallelRead = false;
					if(DigitalNVM*temp = dynamic_cast<DigitalNVM*>(arrayIH->cell[0][0]))
					{    digitalNVM = true;
						if(static_cast<DigitalNVM*>(arrayIH->cell[0][0])->parallelRead == true) parallelRead = true; 
					}

					if(digitalNVM && parallelRead) // parallel read-out for DigitalNVM
					{
						/* TODO: insert algorithm using digitalNVM and parallel readout */
					}
					else {	 // Digital NVM or SRAM row-by-row readout				
						int Dsum = 0;
						int Dsum_logic = 0; 

						for (int k=0; k<param->problemSizeCols; k++) {
							if ((HWInput[k]>>n) & 1) {    // if the nth bit of HWInput[k] is 1
								Dsum += (int)(arrayIH->ReadCell(j,k));
								Dsum_logic += (int)(arrayIH->ReadCell(j,k))*pow(2,n);
							}
						}
						if (DigitalNVM *temp = dynamic_cast<DigitalNVM*>(arrayIH->cell[0][0])) {    // Digital eNVM
							sumArrayReadEnergyIH  += static_cast<DigitalNVM*>(arrayIH->cell[0][0])->readEnergy * arrayIH->numCellPerSynapse * arrayIH->arrayRowSize;
						} 
						else {    // SRAM
							sumArrayReadEnergyIH  += static_cast<SRAM*>(arrayIH->cell[0][0])->readEnergy * arrayIH->numCellPerSynapse * arrayIH->arrayRowSize;
						}

						// Once energy calculation is complete, the value of Dsum we care about for algorithm purposes is stored in Dsum_logic. 
						Dsum = Dsum_logic; 

						// The computed current is the one corresponding to column j, hence to output j.
						Iout[j] = (Dsum);

						// Update maximum column current if needed. 
						if (Iout[j] > maxIout) maxIout = Iout[j];						
					}
				}
			}
		}

		printf("\tMVM Complete. Maximum current: %.4e A\n\tInitiating output extraction.\n",maxIout); 

		/* ********************************************
		 * 				Extraction of Output
		 * 
		 * The output is quantized with numBitPartialSum
		 * precision, and using the maximum current read
		 * as full-scale reference. 
		 * ********************************************/
		if (maxIout == 0) maxIout = 1;
		for (int j = 0; j < param->problemSizeCols; j++) {
			Output[j] = quantize(Iout[j]/maxIout,param->numBitPartialSum)*maxIout; 
	//		printf("\tQuantized Output[%d] = %.4e\n",j,Output[j]);
		}

		// Reset full scale reference for next iteration.
		maxIout = 0; 

		printf("\tOutput normalization phase.\n"); 

		double maxOutput = 0; 
		if (i == 0) {

			for (int j = 0; j < param->problemSizeCols; j++) {
				Output[j] = (Output[j] + 1) / Ith; 

				if (abs(Output[j]) > maxOutput) maxOutput = abs(Output[j]); 
			}

			if (maxOutput > 1) {
				for (int j = 0; j < param->problemSizeCols; j++) {
					Output[j] = Output[j] / maxOutput; 
				}
			}
		}
		else {

			for (int j = 0; j < param->problemSizeCols; j++) {
				Output[j] = Output[j] / Ith; 

				if (abs(Output[j]) > maxOutput) maxOutput = abs(Output[j]);
			}

			if (maxOutput > 1) {
				for (int j = 0; j < param->problemSizeCols; j++) {
					Output[j] = Output[j] / maxOutput; 
				}
			}
		}

		
		/* ********************************************
		 * 					Feedback
		 * 
		 * Output in this time step is re-applied to
		 * the network in the next time step. 
		 * ********************************************/

		for (int j = 0; j < param->problemSizeCols; j++) {
			Input[j] = Output[j];
		}

		/* ********************************************
		 * 					Report
		 * 
		 * After each iteration, we compute the accuracy
		 * of the solution and update the energy report.
		 * ********************************************/

		// Mean absolute error: MAE = sum(|x-x_hat|)/N
		MAE = 0; 
		for (int j = 0; j < param->problemSizeCols; j++) {
			MAE += abs(Output[j] - ev_exact[j]/max_ev);
		}
		MAE = MAE / param->problemSizeCols;

		// Energy report
		numBatchReadSynapse = (int)ceil((double)param->nHide/param->numColMuxed);
		#pragma omp critical    // Use critical here since NeuroSim class functions may update its member variables
		for (int j=0; j<param->nHide; j+=numBatchReadSynapse) {
			int numActiveRows = 0;  // Number of selected rows for NeuroSim
			for (int n=0; n<param->numBitInput; n++) {
				for (int k=0; k<param->nInput; k++) {
					if ((HWInput[k]>>n) & 1) {    // if the nth bit of HWInput on row k is 1
						numActiveRows++;
					}
				}
			}
			subArrayIH->activityRowRead = (double)numActiveRows/param->nInput/param->numBitInput;
			sumNeuroSimReadEnergyIH += NeuroSimSubArrayReadEnergy(subArrayIH);
			sumNeuroSimReadEnergyIH += NeuroSimNeuronReadEnergy(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
			sumReadLatencyIH += NeuroSimSubArrayReadLatency(subArrayIH);
			sumReadLatencyIH += NeuroSimNeuronReadLatency(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
		}

		printf("\tAccuracy at cycle #%d: \t%.2e%\n", i*param->interNumEpochs, MAE);
		printf("\tRead latency: \t%.4e s\n", sumReadLatencyIH);
		printf("\tRead energy: \t%.4e J\n", sumNeuroSimReadEnergyIH);
	}	


	printf("Eigenvector algorithm completed successfully.\n\tv_comp = ["); 
	for(int j = 0; j < param->problemSizeCols; j++) {
		printf("%.4e,",j,Output[j]);
	}
	printf("]';\n");


	printf("\n-----------------------------------------------------------------\n");
	printf("Mean Absolute Error (MAE): \t\t %.4e\n",MAE);
	printf("\nComputation latency, total: \t\t %.4e us\n",sumReadLatencyIH/1e-6);	
	printf("\nComputational performance: \t\t %.4e TOPS\n",pow(param->problemSizeCols,2)/(sumReadLatencyIH)/1e12);	

	printf("\nComputation energy, memory array: \t %.4e pJ\n",sumArrayReadEnergyIH/1e-12);	
	printf("Computation energy, peripherals: \t %.4e pJ\n",sumNeuroSimReadEnergyIH/1e-12);	
	printf("Computation energy, total: \t\t %.4e pJ\n",(sumNeuroSimReadEnergyIH + sumArrayReadEnergyIH)/1e-12);	

	printf("\nEnergy Performance: \t\t\t %.4e TOPS/W\n",pow(param->problemSizeCols,2)/(sumNeuroSimReadEnergyIH + sumArrayReadEnergyIH)/(1e12));	
}
