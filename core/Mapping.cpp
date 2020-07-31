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
#include <string>
#include <vector>
#include <random>
#include "Param.h"
#include "formula.h"
#include "Array.h"
#include "NeuroSim.h"
#include "Debug.h"

extern Param *param;

/* Weights initialization */
void WeightInitialize(std::string mat_name, Array *array) {
    srand(2);

    std::vector<std::vector<double>>    weights(array->arrayColSize, 
                                            std::vector<double>(array->arrayRowSize));


    std::ifstream matrix ("experiments/" + param->experimentName + "/" + mat_name);
    std::string row; 

    /* Prepare counters for maximum and minimum weight */
    param->maxWeight = 0; 
    param->minWeight = INFINITY; 

    /* Initialize weights for the input layer */
    for (int i = 0; i < array->arrayColSize; i++) {
        getline(matrix,row); 
        auto rown = explode(row,';');

        TRACE("\tMatrix initialization: read row %d: [",i);
        for(int jj = 0; jj < rown.size(); jj++) {
            TRACE("%.3f,",stod(rown[jj])); 
        }
        TRACE("]\n");

        for (int j = 0; j < rown.size(); j++) {
            TRACE("%s === %f\n",rown[j].c_str(),stod(rown[j]));

            weights[i][j] = stod(rown[j]);

            param->maxWeight = MAX(param->maxWeight,weights[i][j]);
            param->minWeight = MIN(param->minWeight,weights[i][j]);
        }

        TRACE("\tSuccessfully read matrix file.\n");
        TRACE("\tMinimum weight: %.4e\n",param->minWeight); 
        TRACE("\tMaximum weight: %.4e\n",param->maxWeight); 
        

    }
    matrix.close(); 

    TRACE("\tMatrix initialization complete.\n");

    /* Erase the weight of arrayIH */
    for (int col=0; col<array->arrayColSize; col++) {
        for (int row=0; row<array->arrayRowSize; row++) {
            array->WriteCell(col, row, -(param->maxWeight-param->minWeight), 0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */, param->maxWeight, param->minWeight, false);
        }
    }

    for (int col=0; col<array->arrayColSize; col++) {
        for (int row=0; row<array->arrayRowSize; row++) {   
            array->WriteCell(col, row, weights[col][row], weights[col][row], param->maxWeight, param->minWeight, param->arrayWriteType);
        }
    }
}

/* Conductance initialization (map weight to RRAM conductance or SRAM data) */
void WeightToConductance(Array *array) {
    return;   
}
