//
// Validator.h
// Compares a vessel model with a reconstructed vesselness volume to validate reconstruction.
// 
// Copyright(C) 2024 Sarah F. Frisken, Brigham and Women's Hospital
// 
// This code is free software : you can redistribute it and /or modify it under
// the terms of the GNU General Public License as published by the Free Software 
// Foundation, either version 3 of the License, or (at your option) any later version.
// 
// This code is distributed in the hope that it will be useful, but WITHOUT ANY 
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
// 
// You may have received a copy of the GNU General Public License along with this 
// program. If not, see < http://www.gnu.org/licenses/>.
// 

#pragma once

#include "Volume.h"

class ValidationMetrics
{
public:
	ValidationMetrics() :
		m_TP(0),
		m_TN(0),
		m_FP(0),
		m_FN(0),
		m_totalHaussdorff(0),
		m_numHaussdorffValues(0)
	{};
	~ValidationMetrics() {};

	float truePos() { return m_TP; }
	float trueNeg() { return m_TN; }
	float falsePos() { return m_FP; }
	float falseNeg() { return m_FN; }
	float avgHaussdorff() { return (m_numHaussdorffValues <= 0) ? 0 : m_totalHaussdorff / m_numHaussdorffValues; }
	void incrementTruePos() { m_TP++; }
	void incrementTrueNeg() { m_TN++; }
	void incrementFalsePos() { m_FP++; }
	void incrementFalseNeg() { m_FN++; }
	void incrementHaussdorff(float haussdorfDist) { m_totalHaussdorff += haussdorfDist; m_numHaussdorffValues++; }

private:
	int m_TP;
	int m_TN;
	int m_FP;
	int m_FN;
	float m_totalHaussdorff;
	int m_numHaussdorffValues;
};

class ValidationStats
{
public:
	ValidationStats();
	~ValidationStats();

	float maxDist() { return m_maxDist; }
	float specificity();
	float sensitivity();
	float clDice();
	float accuracy();
	float averageHaussdorff();

	friend class Validator;

private:
	float m_maxDist;
	float m_truePos;
	float m_falsePos;
	float m_trueNeg;
	float m_falseNeg;
	float m_avgHaussdorff;
};
class Validator
{
public:
	Validator(Volume& distVol, Volume& vessselnessVol);
	~Validator();

	ValidationStats validate(float maxDist);

private:
	Volume& m_distVol;
	Volume& m_vesselnessVol;

	void initializeStats(ValidationStats& validationStats, float maxDist);
	void computeStats(ValidationStats& validationStats, float maxDist);
};