//
// Volume.h
// Volume of floating point values
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

#include <array>
#include <QVector3D>

class Volume
{
public:
	Volume(std::array<int, 3> size, float voxelSize, QVector3D offset, float initialValue = 0);
	~Volume();
	void clear(float clearValue = 0);

	// Volume axis is assumed to be aligned with the world coordinate axis. Offset is the position
	// of the volume origin in world coordinates.
	std::array<int,3> size() { return m_size; };
	float voxelSize() { return m_voxelSize; };
	QVector3D origin() { return  m_offset; };

	// Data queries return false if out of bounds
	bool dataAtVoxel(std::array<int, 3> voxel, float& data);	// False if out of bounds
	bool setDataAtVoxel(std::array<int, 3> voxel, float data);	// False if out of bounds
	bool dataAtWorldPoint(QVector3D worlsPoint, float& data);	// Interpolates data at given point. False if out of bounds.
	float* data() { return m_data; };							// Provided for fast data access. No bounds checking.

	// Export volumes as raw NRRD files
	void exportAsScalarNRRD(const char* fileName, bool invert = false);
	void exportAsScalarNRRD(const char* fileName, float minValue, float maxValue, bool invert = false);

private:
	std::array<int, 3> m_size;
	float m_voxelSize;
	QVector3D m_offset;
	float* m_data;

	void writeHeaderFile(const char* filename, std::array<int, 3> size, QVector3D axis[3], QVector3D origin);
	void writeDataFile(const char* filename, std::array<int, 3> size, const char* data);
};