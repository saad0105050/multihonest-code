#pragma once
#include "ReachAndMargin.h"

class ReachAndMarginMultiHonest : public ReachAndMargin
{
	double *prUniqueH, *prMultiH;
	
	void probMassIncomingToReachZeroMarginZero();
	void probMassIncomingToReachZeroMarginNegOne();
	virtual void calcNewValue(const int rho, const int mu);
public:
	ReachAndMarginMultiHonest(
		const int N,
		const int R,
		const DoubleVector& advFractions,
		const double fracUniqueHonest = 1.0
		);
	~ReachAndMarginMultiHonest();
private:
};