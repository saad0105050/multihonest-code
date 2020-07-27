#include "MultiHonest.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <chrono>

struct Params {
	// you can assign the members here
	bool useLogarithms = false;
	bool useTimer = false;
	bool printScientificFormat = false;
	//
	//--------------------------------
	// DO NOT manually assign the members below
	bool everythingOkay = false;
	int N;
	int R;
	DoubleVector advFracArr;
	UintVector Ns;
	double fracUniqueH = 1.0; // pr[uniquely honest slot] = (1 - alpha) * fracUniqueH

	~Params(){
		advFracArr.clear();
		Ns.clear();
	}
};



void write(ReachAndMargin* const rm, const int t, std::ofstream& file, 
	const bool useLogarithms = false, 
	const bool useScientific = false) {
	DoubleVector prForkable = rm->forkableProbability();
	DoubleVector logPrForkable;
	std::stringstream ss;

	// gather values
	if (useLogarithms) {
		for (auto i = 0u; i < prForkable.size(); i++)
			logPrForkable.push_back(std::log10(prForkable[i]));
	}
	DoubleVector& probs = useLogarithms ? logPrForkable : prForkable;
	ss << t << ",\t";
	if (useScientific){
		for (auto p : probs) {
			char dest[100];
			// sprintf_s(dest, "%1.2E, ", p);
			snprintf(dest, 100, "%1.2E, ", p);
			ss << dest;
		}
	}
	else
		ss << probs;
		//<< " (total: " << rm->totalProbability() << ") "
	ss << std::endl;
	//std::cout << ss.str(); 
	file << ss.str(); ss.str("");

}


Params* init(std::vector<double> advFracArr, std::vector<int> NminMaxInc, int R, std::ofstream& file) {

	Params* pParams = new Params;
	Params& params = *pParams;
	std::stringstream ss;

	//const double advFracArr[] = { 0.05, 0.1};
	params.advFracArr = advFracArr;
	if (params.advFracArr.size() > ReachAndMargin::MAX_FRACS) {
		ss << "Maximum number of adversarial fractions: " << ReachAndMargin::MAX_FRACS
			<< " but we have: " << params.advFracArr.size() << std::endl;
		std::cerr << ss.str(); file << ss.str(); ss.str("");
		params.everythingOkay = false;
		return pParams;
	}

	auto Nmin = NminMaxInc[0],
		Nmax = NminMaxInc[1],
		Ninc = NminMaxInc[2];
	for (int n = Nmin; n <= Nmax; n += Ninc)
		params.Ns.push_back(n);
	params.N = *std::max_element(params.Ns.begin(), params.Ns.end());
	if (params.N > ReachAndMargin::NMAX) {
		ss << "Maximum possible N: " << ReachAndMargin::NMAX
			<< " but we have: " << params.N << std::endl;
		std::cerr << ss.str(); file << ss.str(); ss.str("");
		params.everythingOkay = false;
		return pParams;
	}

	params.R = R;
	if (R < 0 || R > ReachAndMargin::RMAX) {
		ss << "R must be between 0 and " << ReachAndMargin::RMAX << " (inclusive)" 
			<< " but we have: " << R << std::endl;
		std::cerr << ss.str(); file << ss.str(); ss.str("");
		params.everythingOkay = false;
		return pParams;
	}
	params.everythingOkay = true;
	return pParams;
}


void print_header(Params& params, ReachAndMargin* const rm, std::ofstream& file){
	std::stringstream ss;
	// Print header
	ss << "\n\nCreating instance with N = " << params.N << ", R = " << params.R
		<< ", and fracUniqueH: " << params.fracUniqueH 
		<< ", and advFracs: " << params.advFracArr << std::endl;

	ss << "Using memory: " << rm->totalMemory() << " bytes" << std::endl;
	ss << std::endl << "N," << "\t" << params.advFracArr << std::endl;
	std::cout << ss.str();
	file << ss.str();
	// clean up
	ss.clear();
}

void run_experiment(Params& params, std::ofstream& file){
	if (!params.everythingOkay) return;

	// create main object
	auto rm = new ReachAndMarginMultiHonest(params.N, params.R, params.advFracArr, params.fracUniqueH);
	params.R = rm->R;

	std::cout << "R: " << params.R;
	std::cout << std::endl;

	// print header in file
	print_header(params, rm, file);

	int toBeWritten = 0;

	std::chrono::system_clock::time_point start;
	std::chrono::duration<double> sum, last;

	// evolve state
	for (int t = 1; t <= params.N; t++) {
		std::cout << "\rEvolution: " << std::setw(4) << t << " of " << params.N;
		if (params.useTimer) {
			std::cout 
				<< "; last evo: " << std::setw(10) << last.count() << " sec;"
				<< " total: " << std::setw(10) << sum.count() << " sec\n" << std::flush;
		}

		if (params.useTimer) {
			start = std::chrono::high_resolution_clock::now();
		}
		rm->evolve();

		if (params.useTimer) {
			auto finish = std::chrono::high_resolution_clock::now();
			last = finish - start;
			sum += last;
		}
		// write current state to file
		if (t == (int) params.Ns[toBeWritten]) {
			write(rm, t, file, 
				params.useLogarithms, 
				params.printScientificFormat);
			toBeWritten++;
		}
	}
	// Clean up
	DELETE(rm);
}

void make_forkability_table(
	std::vector<double> advFracArr, 
	std::vector<int> NminMaxInc, 
	std::vector<double> fracUniqueHonestArr,
	int R, std::string fileName, 
	const bool useLogarithm = true,
	const bool printScientificFormat = false
	){

	// Test whether memory leak detection is working
	//auto test_mem_leak = new std::string("Testing mem leak detection"); // this will cause memory leaks unless deleted

#ifdef _WIN32   
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	std::ofstream file;
	file.open(fileName);


	auto pParams = init(advFracArr, NminMaxInc, R, file);
	pParams->useLogarithms = useLogarithm;
	pParams->printScientificFormat = printScientificFormat;

	// one run for each prUniqueHonest
	for (auto f : fracUniqueHonestArr) {
		pParams->fracUniqueH = f;
		if (pParams->everythingOkay) 
			run_experiment(*pParams, file);
	}

	// finally
	file.close();

	// display on terminal
	std::cout << "\n------------- Output (" << fileName << ") ------------" << std::endl;
	std::stringstream ss;
	ss << "cat " << fileName;
	system(ss.str().c_str());

	// clean up
	DELETE(pParams);
	file.clear();
	ss.clear();
	fileName.clear();


#ifdef _WIN32   
	_CrtDumpMemoryLeaks();
#endif

}

// See the definition of struct Params for options you can set/unset
int main() {
	const char* default_fileName = "forkability_probs.multihonest.txt";
	std::string fileName;

	std::cout << "Enter output filename [" << default_fileName << "]: ";
	std::getline(std::cin, fileName);
	if (fileName.empty())
		fileName = std::string(default_fileName);
	std::cout << "Using output file: " << fileName << std::endl;

	// Adversarial fractions
	std::vector<double> advFracArr = { 
		// 0.01, 0.49 
		0.01, 0.1, 0.2, 0.3, 0.4, 0.49
	};

	// Array {MIN, MAX, INCREMENT}
	// Timeslots at which to write Pr[forkable] in file,
	// starting from MIN to MAX with INCREMENT
	// Note that the Markov chain will be evolved up to MAX timeslots
	std::vector<int> NminMaxInc = { 
		// 10, 10, 10
		100, 500, 100
	};

	// Vector of ph/(1 - pA) values
	std::vector<double> fracUniqueHonestArr = { 
		// 0.01, 0.25
		0.01, 0.25, 0.5, 0.8, 0.9, 1.0 
	};


	// Support of the initial reach distribution will be [0, min(R, RMAX)]
	int R = ReachAndMargin::RMAX;  
	//int R = 10;

	bool	useLogarithm = false,
			printScientificFormat = true
			;



	// do the work
	make_forkability_table(
		advFracArr, NminMaxInc, fracUniqueHonestArr,
		R, fileName, 
		useLogarithm, printScientificFormat
	);


	fileName.clear();
	// system("pause");
	return 0;
}