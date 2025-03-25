//
// BravaUtils.h
// Tools for reading and parsing vessel data from Brava .swc files.
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

#include <vector>

#include "VesselModel.h"

class BravaUtils
{
public:
	BravaUtils();
	~BravaUtils();

	void importBravaData(const char* filename);

	// Get vessel vertices for the specified circulation. Brava circulations are 
	// encoded as: 0 = COW, 1 = none, 2 = RightPosteriorCA, 3 = RightMiddleCA, 
	// 4 = LeftMiddleCA, 5 = LeftPosteriorCA, 6 = RightAnteriorCA, 7 = LeftAnteriorCA
	// (COW is Circle of Willis, CA is cerebral artery)
	std::vector<VesselModel::Line> vesselLines(std::vector<int> circulations);

private:
	struct Branch {
		int ID;
		std::vector<VesselModel::Vertex> vertices;
	};
	std::vector<Branch> m_branches;
};