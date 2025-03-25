//
// Reconstructor.h
// Reconstructs a 3D vesselness volume from a set of annotated 2D projection images.
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

#include "ProjectionImage.h"
#include "Validator.h"
#include "Volume.h"

#include <string>
#include <vector>

class ReconstructionThresholds
{
public:
	ReconstructionThresholds(float maxDistThreshold = 0, float maxRadiusDiffThreshold = 0, 
		float maxTimeDiffThreshold = 0, float maxSharedVelCompDiffThreshold = 0) :
		m_tDist(maxDistThreshold),
		m_tRadius(maxRadiusDiffThreshold),
		m_tTime(maxTimeDiffThreshold),
		m_tVelocityComponent(maxSharedVelCompDiffThreshold)
	{
	}
	~ReconstructionThresholds() {};

	float maxDist() { return m_tDist; };
	float maxRadiusDiff() { return m_tRadius; };
	float maxTimeDiff() { return m_tTime; };
	float maxSharedVelCompDiff() { return m_tVelocityComponent; };

private:
	float m_tDist;
	float m_tRadius;
	float m_tTime;
	float m_tVelocityComponent;
};

class ReconstructionParameters
{
public:
	enum class ReconstructionType { MaxDist, AvgDist, ConstrainedDist };

	ReconstructionParameters(ReconstructionType type = ReconstructionType::MaxDist, float maxDist = 10, 
		float voxelSize = 1, std::array<int, 3> volumeSize = { 0, 0, 0 }, 
		QVector3D volumeOrigin = { 0, 0, 0 }) :
		m_type(type),
		m_maxDist(maxDist),
		m_voxelSize(voxelSize),
		m_volumeSize(volumeSize),
		m_volumeOrigin(volumeOrigin)
	{
	};
	~ReconstructionParameters() {};

	ReconstructionType type() { return m_type; }
	float maxDist() { return m_maxDist; }
	float voxelSize() { return m_voxelSize; }
	std::array<int, 3> volumeSize() { return m_volumeSize; }
	QVector3D volumeOrigin() { return m_volumeOrigin; }

private:
	ReconstructionType m_type;
	float m_maxDist;
	float m_voxelSize;
	std::array<int, 3> m_volumeSize;	// Dimensions
	QVector3D m_volumeOrigin;			// Left, lower, back corner of axis-aligned volume in world coordinates
};

class Reconstructor
{
public:
	Reconstructor(std::vector<std::shared_ptr<ProjectionImage>> projImages, ReconstructionParameters reconParams);
	~Reconstructor();

	void reconstructAndExport(const char* filename);
	void reconstructAndExport(const char* filename, ReconstructionThresholds thresholds);
	std::vector<ValidationMetrics> reconstructAndValidate(std::shared_ptr<Volume> distField,
		std::vector<ReconstructionThresholds> thresholds);
	ValidationMetrics reconstructAndValidateMaxDist(std::shared_ptr<Volume> distField, float maxDistThreshold);

private:
	std::vector<std::shared_ptr<ProjectionImage>> m_projImages;
	ReconstructionParameters m_reconParams;

	void reconstructAndExportAvgDist(const char* filename);
	void reconstructAndExportMaxDist(const char* filename);
	void reconstructAndExportConstrainedDist(const char* filename, ReconstructionThresholds thresholds);
	std::vector<ValidationMetrics> reconstructAndValidateConstrainedDist(std::shared_ptr<Volume> distField,
		std::vector<ReconstructionThresholds> thresholds);

	void testProjection();
};