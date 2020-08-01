/*******************************************************************************
* Copyright (c) 2020
* Dipartimento di Elettronica, Informazione e Bioingegneria, Politecnico di Milano
* PI: Prof. Daniele Ielmini
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of SpikUSim - a simulator to benchmark the SpikUS 
* architecture based on NeuroSim, a device-circuit-algorithm framework to benchmark 
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
*   Piergiulio Mannocci     Email: piergiulio dot mannocci at polimi dot it 
*
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
*                     
********************************************************************************/

#include <string>
#include "math.h"
#include "Results.h"
#include "Param.h"
#include "formula.h"

/* Logger libraries */
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

extern spdlog::logger logger;

extern Param* param; 

Results::Results() {
	
	MAE = INFINITY; 
	MRE = INFINITY; 

    throughput = 0;  
	energyPerformance = 0; 
	performanceDensity = 0; 
	SWaP = 0; 
	
	totalEnergy = 0; 
	totalCoreEnergy = 0; 
	totalNeuronEnergy = 0;
	totalCorePeripheralEnergy = 0; 
	totalCoreMemoryEnergy = 0; 

	totalLatency = 0; 

	totalArea = 0; 
	totalCoreArea = 0; 
	totalNeuronArea = 0; 
	totalCorePeripheralArea = 0; 
	totalCoreMemoryArea = 0;

	totalNeuronLeakage = 0; 
	totalCoreLeakage = 0; 
	totalFlops = 0; 

}

/** Compute benchmarks for the currently stored 
 * values of energy, latency and area.
 * In case an equivalent benchmark is sought, 
 * computed values are adapted following the 
 * adaption parameters. 
 * */
void Results::Benchmark() {

	/* Adapt stored energy, latency and area values if needed. */
	totalCoreMemoryArea *= GetAdapter(Category::AREA,Element::MEMORY);
	totalCorePeripheralArea *= GetAdapter(Category::AREA,Element::PERIPHERALS);
	totalNeuronArea *= GetAdapter(Category::AREA,Element::NEURONS);

	totalCoreMemoryEnergy *= GetAdapter(Category::ENERGY,Element::MEMORY); 
	totalCorePeripheralEnergy *= GetAdapter(Category::ENERGY,Element::PERIPHERALS);
	totalNeuronEnergy *= GetAdapter(Category::ENERGY,Element::NEURONS);

	totalLatency *= GetAdapter(Category::LATENCY,Element::SYSTEM);


	totalCoreArea = totalCoreMemoryArea + totalCorePeripheralArea; 
	totalArea = totalCoreArea + totalNeuronArea; 

	totalCoreEnergy = totalCoreMemoryEnergy + totalCorePeripheralEnergy;
	totalEnergy = totalCoreEnergy + totalNeuronEnergy;

	totalFlops = pow(param->problemSize,2)*param->numCycles;
	throughput = totalFlops/totalLatency;
	energyPerformance = totalFlops/totalEnergy;
	performanceDensity = throughput/totalArea; 
	SWaP = energyPerformance/totalArea; 
}

/** Generate CSV file containing benchmark results. */
void Results::GenerateCSV() {

}

/** Print benchmark results to the console. 
 * */
void Results::PrintResults() {

	std::cout << "Benchmark results:" << std::endl; 
	std::cout << "MAE: " << MAE << std::endl; 
	std::cout << "MRE: " << MRE << std::endl; 
	std::cout << "Energy: " << totalEnergy/1e-9 << " nJ" << std::endl; 
	std::cout << "Latency: " << totalLatency/1e-6 << " us" << std::endl; 
	std::cout << "Area: " << totalArea*1e6 << " mm^2" << std::endl; 
	std::cout << std::endl; 
	std::cout << "Throughput: " << throughput/1e9 << " GOPS" << std::endl; 
	std::cout << "Performance: " << energyPerformance/1e12 << " TOPS/W" << std::endl; 
	std::cout << "Density: " << performanceDensity/(1e12*1e6) << " TOPS/mm^2" << std::endl; 
	std::cout << "SWaP: " << SWaP/(1e12*1e6) << " TOPS/(W mm^2)" << std::endl; 

	std::cout << std::endl << std::endl; 
	std::cout << "Coarse breakdown:" << std::endl; 
	std::cout << "Memory energy: " << totalCoreMemoryEnergy/1e-9 << " nJ" << std::endl; 
	std::cout << "Peripherals energy: " << totalCorePeripheralEnergy/1e-9 << " nJ" << std::endl; 
	std::cout << "Neuron energy: " << totalNeuronEnergy/1e-9 << " nJ" << std::endl; 
	std::cout << std::endl; 
	std::cout << "Memory area: " << totalCoreMemoryArea*1e6 << " mm^2" << std::endl; 
	std::cout << "Peripherals area: " << totalCorePeripheralArea*1e6 << " mm^2" << std::endl; 
	std::cout << "Neuron area: " << totalNeuronArea*1e6 << " mm^2" << std::endl; 

	return; 
}

/** Print benchmark results to an output file. 
 * @param outstr 	(const char*) the output file.
 * */
void Results::PrintResultsToFile(const char* outstr) {
	std::ofstream outfile(outstr);

	outfile << "Benchmark results:" << std::endl; 
	outfile << "MAE: " << MAE << std::endl; 
	outfile << "MRE: " << MRE << std::endl; 
	outfile << "Energy: " << totalEnergy/1e-9 << " nJ" << std::endl; 
	outfile << "Latency: " << totalLatency/1e-6 << " us" << std::endl; 
	outfile << "Area: " << totalArea*1e6 << " mm^2" << std::endl; 
	outfile << std::endl; 
	outfile << "Throughput: " << throughput/1e9 << " GOPS" << std::endl; 
	outfile << "Performance: " << energyPerformance/1e12 << " TOPS/W" << std::endl; 
	outfile << "Density: " << performanceDensity/(1e12*1e6) << " TOPS/mm^2" << std::endl; 
	outfile << "SWaP: " << SWaP/(1e12*1e6) << " TOPS/(W mm^2)" << std::endl; 

	outfile << std::endl << std::endl; 
	outfile << "Coarse breakdown:" << std::endl; 
	outfile << "Memory energy: " << totalCoreMemoryEnergy/1e-9 << " nJ" << std::endl; 
	outfile << "Peripherals energy: " << totalCorePeripheralEnergy/1e-9 << " nJ" << std::endl; 
	outfile << "Neuron energy: " << totalNeuronEnergy/1e-9 << " nJ" << std::endl; 
	outfile << std::endl; 
	outfile << "Memory area: " << totalCoreMemoryArea*1e6 << " mm^2" << std::endl; 
	outfile << "Peripherals area: " << totalCorePeripheralArea*1e6 << " mm^2" << std::endl; 
	outfile << "Neuron area: " << totalNeuronArea*1e6 << " mm^2" << std::endl; 

	outfile.close(); 
	return; 
}

double Results::GetAdapter(Category category, Element element) {

	auto ads = explode(param->adaptBenchmark,'|');

	for (int i = 0; i < ads.size(); i++) {

		auto cat = explode(ads[i],':');

		if (cat[0].c_str()[0] != ((char)category)) continue; 

		auto el_ad_pairs = explode(cat[1],';');

		for (int j = 0; j < el_ad_pairs.size(); j++) {

			auto el = explode(el_ad_pairs[j],'*');

			if (el[0].c_str()[0] != ( (char) element) ) continue; 

			logger.warn("Adapted {},{} by {}",category,element,el[1]);
			return stod(el[1]);
		}
	}

	return 1;
}