#include "MultiHonest.h"
#include<iostream>

ReachAndMarginMultiHonest::ReachAndMarginMultiHonest(
		const int N, const int R, 
		const DoubleVector& advFractions,
		const double fracUniqueHonest)
		: ReachAndMargin(N, R, advFractions) {
	prUniqueH = new double[NFRACS];
	prMultiH = new double[NFRACS];
	for (int i = 0; i < NFRACS; i++) {
		prUniqueH[i] = prDown[i] * fracUniqueHonest;
		prMultiH[i] = prDown[i] - prUniqueH[i];
	}

}

ReachAndMarginMultiHonest::~ReachAndMarginMultiHonest() {
	DELETE_ARR(prUniqueH);
	DELETE_ARR(prMultiH);
}

void ReachAndMarginMultiHonest::probMassIncomingToReachZeroMarginZero() {
	const int r = 0, m = 0;

	// from (r+1, m) using downstep, 
	// from (r+1, m+1) using downstep
	auto right = get(prevMat, time - 1, r + 1, m);
	auto diagonallyUp = get(prevMat, time - 1, r + 1, m + 1);

	add(tempA, diagonallyUp, right);
	mult(tempB, prDown, tempA);

	// from (r, m) using multihonest downstep
	auto same_cell = get(prevMat, time - 1, r, m);
	mult(tempA, prMultiH, same_cell);
	// finally
	auto newValue = get(curMat, time, r, m);
	add(newValue, tempA, tempB);
}

void ReachAndMarginMultiHonest::probMassIncomingToReachZeroMarginNegOne() {
	const int r = 0, m = -1;
	// from (r, m+1) using uniquely honest downstep
	auto up = get(prevMat, time - 1, r, m + 1);
	auto newValue = get(curMat, time, r, m);
	mult(newValue, prUniqueH, up);
}

void ReachAndMarginMultiHonest::calcNewValue(const int r, const int m) {
	// initialize to zero
	//for (int a = 0; a < NFRACS; a++)
	//	newValue[a] = 0;

	if (r >= 1 || (m <= -2 || m >= 1) )
		// not a multihonest case
		ReachAndMargin::calcNewValue(r, m);
	else if (r == 0 && m == 0) {
		// this is a special multihonest case
		// margin stays the same with prMultiH, 
		probMassIncomingToReachZeroMarginZero();
	}
	else if (r == 0 && m == -1) {
		// this is a special multihonest case
		// it can be reached from (r = 0, m = 0) with prUniqueH
		probMassIncomingToReachZeroMarginNegOne();
	}
}
