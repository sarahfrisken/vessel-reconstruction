//
// Validator.cpp
// Implementation of Validator.
//

#include "Validator.h"

// 
// Public
//
ValidationStats::ValidationStats() :
	m_maxDist(0),
	m_truePos(0),
	m_falsePos(0),
	m_trueNeg(0),
	m_falseNeg(0),
	m_avgHaussdorff(0)
{
}
ValidationStats::~ValidationStats()
{
}

float ValidationStats::specificity()
{
	float specificity = 0;
	if (float denominator = m_trueNeg + m_falsePos > 0) {
		specificity = m_trueNeg / denominator;
	}
	return specificity;
}
float ValidationStats::sensitivity()
{
	float sensitivity = 0;
	if (float denominator = m_truePos + m_trueNeg > 0) {
		sensitivity = m_truePos / denominator;
	}
	return sensitivity;
}
float ValidationStats::clDice()
{
	float clDice = 0;
	if (float denominator = 2 * m_truePos + m_falsePos + m_falseNeg > 0) {
		clDice = 2 * m_truePos / denominator;
	}
	return clDice;
}
float ValidationStats::accuracy()
{
	float accuracy = 0;
	if (float denominator = m_truePos + m_trueNeg + m_falsePos + m_falseNeg > 0) {
		accuracy = (m_truePos + m_trueNeg) / denominator;
	}
	return accuracy;
}
float ValidationStats::averageHaussdorff()
{
	return m_avgHaussdorff;
}

// 
// Public
//
Validator::Validator(Volume& distVol, Volume& vessselnessVol) :
	m_distVol(distVol),
	m_vesselnessVol(vessselnessVol)
{

}
Validator::~Validator()
{
}

ValidationStats Validator::validate(float maxDist)
{
	ValidationStats validationStats;
	initializeStats(validationStats, maxDist);
	computeStats(validationStats, maxDist);
	return validationStats;
}

//
// Private
//
void Validator::initializeStats(ValidationStats& validationStats, float maxDist)
{
	validationStats.m_maxDist = maxDist;
	validationStats.m_truePos = 0;
	validationStats.m_falsePos = 0;
	validationStats.m_trueNeg = 0;
	validationStats.m_falseNeg = 0;
	validationStats.m_avgHaussdorff = 0;
}
void Validator::computeStats(ValidationStats& validationStats, float maxDist)
{

}


