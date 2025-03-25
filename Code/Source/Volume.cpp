//
// Volume.cpp
// Implementation of Volume.
//

#include "Volume.h"

#include <iostream>
#include <fstream>
// 
// Public
//
Volume::Volume(std::array<int, 3> size, float voxelSize, QVector3D offset, float initialValue) :
	m_size(size),
	m_voxelSize(voxelSize),
	m_offset(offset),
	m_data(nullptr)
{
	size_t arraySize = (size_t) m_size[0] * m_size[1] * m_size[2];
	try {
		m_data = new float[arraySize];
	}
	catch (const std::bad_alloc& e) {
		std::cout << "Distance volume allocation failed: " << e.what() << '\n';
	}

	float* pData = m_data;
	for (int i = 0; i < arraySize; i++) {
		*pData++ = initialValue;
	}
}
Volume::~Volume()
{
	delete[] m_data;
	m_data = nullptr;
}
void Volume::clear(float value)
{
	if (!m_data) return;
	float* pData = m_data;
	int numVoxels = m_size[0] * m_size[1] * m_size[2];
	for (int i = 0; i < numVoxels; i++) *pData++ = value;
}


bool Volume::dataAtVoxel(std::array<int, 3> voxel, float& data)
{
	if (!m_data) return false;
	if (voxel[0] < 0 || voxel[0] >= m_size[0] ||
		voxel[1] < 0 || voxel[1] >= m_size[1] ||
		voxel[2] < 0 || voxel[2] >= m_size[2]) {
		return false;
	}
	data = m_data[voxel[0] + voxel[1] * m_size[0] + voxel[2] * m_size[0] * m_size[1]];
	return true;
}
bool Volume::setDataAtVoxel(std::array<int, 3> voxel, float data)
{
	if (!m_data) return false;
	if (voxel[0] < 0 || voxel[0] >= m_size[0] ||
		voxel[1] < 0 || voxel[1] >= m_size[1] ||
		voxel[2] < 0 || voxel[2] >= m_size[2]) {
		return false;
	}
	m_data[voxel[0] + voxel[1] * m_size[0] + voxel[2] * m_size[0] * m_size[1]] = data;
	return true;
}
bool Volume::dataAtWorldPoint(QVector3D worldPoint, float& data)
{
	if (!m_data || m_voxelSize == 0) return false;
	QVector3D point = (worldPoint - m_offset) / m_voxelSize;
	if (point[0] < 0 || point[0] >= m_size[0] ||
		point[1] < 0 || point[1] >= m_size[1] ||
		point[2] < 0 || point[2] >= m_size[2]) {
		return false;
	}

	// Compute bilinear interpolation parameters
	int i0 = (int)point.x();
	int j0 = (int)point.y();
	int k0 = (int)point.z();
	float dx = point.x() - i0;
	float dy = point.y() - j0;
	float dz = point.z() - k0;
	int i1 = i0 + 1;
	int j1 = j0 + 1;
	int k1 = k0 + 1;
	float ddx = 1 - dx;
	float ddy = 1 - dy;
	float ddz = 1 - dz;

	// Avoid sampling outside the volume
	if (i0 == m_size[0] - 1) { i0 = m_size[0] - 2; i1 = i0 + 1; dx = 0; ddx = 1; };
	if (j0 == m_size[1] - 1) { j0 = m_size[1] - 2; j1 = j0 + 1; dy = 0; ddy = 1; };
	if (k0 == m_size[2] - 1) { k0 = m_size[2] - 2; k1 = k0 + 1; dz = 0; ddz = 1; };

	// Perform binlinear interpolation
	float d000 = m_data[i0 + j0 * m_size[0] + k0 * m_size[0] * m_size[1]];
	float d100 = m_data[i1 + j0 * m_size[0] + k0 * m_size[0] * m_size[1]];
	float d010 = m_data[i0 + j1 * m_size[0] + k0 * m_size[0] * m_size[1]];
	float d110 = m_data[i1 + j1 * m_size[0] + k0 * m_size[0] * m_size[1]];
	float d001 = m_data[i0 + j0 * m_size[0] + k1 * m_size[0] * m_size[1]];
	float d101 = m_data[i1 + j0 * m_size[0] + k1 * m_size[0] * m_size[1]];
	float d011 = m_data[i0 + j1 * m_size[0] + k1 * m_size[0] * m_size[1]];
	float d111 = m_data[i1 + j1 * m_size[0] + k1 * m_size[0] * m_size[1]];

	data = ddx * ddy * ddz * d000 +
		dx  * ddy * ddz * d100 +
		ddx * dy  * ddz * d010 +
		dx  * dy  * ddz * d110 +
		ddx * ddy * dz  * d001 +
		dx  * ddy * dz  * d101 +
		ddx * dy  * dz  * d011 +
		dx  * dy  * dz  * d111;
	return true;
}

void Volume::exportAsScalarNRRD(const char* fileName, bool invert)
{
	if (m_data == nullptr) return;

	// Determine the volume intensity range automatically 
	float minVal = m_data[0];
	float maxVal = m_data[0];
	float* pData = m_data;
	for (int k = 0; k < m_size[2]; k++) {
		for (int j = 0; j < m_size[1]; j++) {
			for (int i = 0; i < m_size[0]; i++) {
				minVal = std::min(minVal, *pData);
				maxVal = std::max(maxVal, *pData);
				pData++;
			}
		}
	}
	exportAsScalarNRRD(fileName, minVal, maxVal, invert);
}
void Volume::exportAsScalarNRRD(const char* fileName, float minValue, float maxValue, bool invert)
{
	float valDif = maxValue - minValue;
	if (m_data == nullptr || valDif <= 0) return;

	// Generate a scaled char* array of intensity values
	char* data = new char[(size_t)m_size[0] * m_size[1] * m_size[2]];
	char* pCharData = data;
	float offsetToChar = minValue;
	float scaleToChar = 255.0 / valDif;
	float* pValue = m_data;
	for (int k = 0; k < m_size[2]; k++) {
		for (int j = 0; j < m_size[1]; j++) {
			for (int i = 0; i < m_size[0]; i++) {
				float scaledToChar = std::min(std::max(0.0f, (*pValue++ - offsetToChar)* scaleToChar), 255.0f);
				char scaledChar = invert ? 255 - (char)scaledToChar : (char)scaledToChar;
				*pCharData++ = scaledChar;
			}
		}
	}

	QVector3D axis[3] = {
		{ m_voxelSize, 0, 0 },
		{ 0, m_voxelSize, 0 },
		{ 0, 0, m_voxelSize }
	};
	if (data == nullptr) return;
	writeHeaderFile(fileName, m_size, axis, m_offset);
	writeDataFile(fileName, m_size, reinterpret_cast<const char*>(data));
}

// 
// Private
//
// Export Raw NRRD file
// Create and write the NRRD header file 
// In the RAS coordinate system, origin is the left, posterior, inferior corner. Axes 
// go from left to right, posterior to anterior, and inferior to superior.
//
// Space directions are input as volume coordinate axes in world space RAS coordinates, 
// with axis[0] corresponding to the image x axis, axis[0] corresponding to the image y 
// axis and axis[2] corresponding to the z axis for volumes or the image normal direction.
// For example, an image in the coronal plane (RS) with origin at the origin is left, 
// posterior, inferior corner would have axes: (1, 0, 0) (0, 0, 1) (0, 1, 0) and origin 
// (0,0,0). The same image offset in the anterior direction by 100 voxels, would have
// the same axes and origin (0, 100, 0). The image rotated by A degrees about the S
// axis with origin (l, a, s) would have axes:  (cosA,sinA,0) (0,0,1) (-sinA,cosA,0)
// and origin (l,a,s).
void Volume::writeHeaderFile(const char* filename, std::array<int, 3> size, QVector3D axis[3],
	QVector3D origin)
{
	std::string hdrFilename = filename;
	hdrFilename.append(".nhdr");
	std::ofstream headerFile(hdrFilename.c_str());
	headerFile << "NRRD0004\n";
	headerFile << "# Complete NRRD file format specification at: \n";
	headerFile << "# http://teem.sourceforge.net/nrrd/format.html \n";
	headerFile << "type: unsigned char\n";
	headerFile << "dimension: 3 \n";
	headerFile << "space: right-anterior-superior\n";
	headerFile << "sizes: " << std::to_string((int)size[0]) << " " << std::to_string((int)size[1]) << " " << std::to_string((int)size[2]) << "\n";
	std::string sAxis0 = " (" + std::to_string(axis[0][0]) + "," + std::to_string(axis[0][1]) + "," + std::to_string(axis[0][2]) + ")";
	std::string sAxis1 = " (" + std::to_string(axis[1][0]) + "," + std::to_string(axis[1][1]) + "," + std::to_string(axis[1][2]) + ")";
	std::string sAxis2 = " (" + std::to_string(axis[2][0]) + "," + std::to_string(axis[2][1]) + "," + std::to_string(axis[2][2]) + ")";
	headerFile << "space directions:" << sAxis0 << sAxis1 << sAxis2 << "\n";
	headerFile << "kinds: space space space\n";
	headerFile << "encoding: raw\n";
	std::string sOrigin = " (" + std::to_string(origin[0]) + "," + std::to_string(origin[1]) + "," + std::to_string(origin[2]) + ")";
	headerFile << "space origin: " << sOrigin << "\n";
	headerFile << "data file: " << filename << ".raw\n";
	headerFile.close();
}

// Create and write the data file 
void Volume::writeDataFile(const char* filename, std::array<int, 3> size, const char* data)
{
	std::string dataFilename = filename;
	dataFilename.append(".raw");
	std::ofstream dataFile(dataFilename, std::ios::binary);
	int dataSize = size[0] * size[1] * size[2];
	dataFile.write(data, dataSize);
	dataFile.close();
}