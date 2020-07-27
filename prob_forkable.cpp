#include "MultiHonest.h"
#include <iostream>

void prob_forkable(const int N, const int R, double advFrac, double fracUniqueH) {
	auto advFracArr = { advFrac };
	//ReachAndMargin* rm = new ReachAndMargin(N, R, advFracArr);
	auto rm = new ReachAndMarginMultiHonest(N, R, advFracArr, fracUniqueH);
	for (int t = 1; t <= N; t++) {
		rm->evolve();
	}
	std::cout << "\tPr[forkable] <= " << rm->forkableProbability() << std::endl;
	DELETE(rm);

}


int main(){
	int N = 10; 
	// int	R = N;	
	double	eps = 0.5,
			fracUniqueH = 1.0;

	std::cout <<
		"N is the length of execution (integer, at least 1, default 10)" << std::endl <<
		// "R is the maximum reach of the prefix (between 0 and " << ReachAndMargin::RMAX << " inclusive, default 10)" << std::endl <<
		"eps is the bias (at least 0 but less than 1, default 0.5)" << std::endl <<
		"fracUniqueH is the ratio of #uniquely honest slots to #honest slots (at least 0 but less than 1, default 1.0)" << std::endl << std::endl;

	std::cout << "Enter N" << 
		// ", R" << 
		", eps, and fracUniqueH (separate by space): "; 
	if (std::cin.peek() != '\n') {
		std::cin >> N;
		// std::cin >> R;
		std::cin >> eps;
		std::cin >> fracUniqueH;
	}
	if (N < 1) {
		std::cerr << "N must be at least 1";
		exit(-1);
	}
	// if (R < 0 || R > ReachAndMargin::RMAX) {
	// 	std::cerr << "R must be between 0 and " << ReachAndMargin::RMAX << "  (inclusive)";
	// 	exit(-1);
	// }
	if (eps < 0 || eps > 1) {
		std::cerr << "eps must be between 0 and 1";
		exit(-1);
	}
	if (fracUniqueH < 0 || fracUniqueH > 1) {
		std::cerr << "fracUniqueH must be between 0 and 1";
		exit(-1);
	}
	std::cout
		<< "Received N: " << N
		// << " R: " << R
		<< " eps: " << eps
		<< " fracUniqueH: " << fracUniqueH
		<< std::endl;
	
	double advFrac = (1 - eps) / 2.0;
	auto R = ReachAndMargin::RMAX;
	prob_forkable(N, R, advFrac, fracUniqueH);

	// system("pause");
	return 0;
}