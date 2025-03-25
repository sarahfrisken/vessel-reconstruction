//
// ProjectionImage.cpp
// Implementation of ProjectionImage.
//

#include "ProjectionImage.h"

#include <iostream>
#include <fstream>

// 
// Public
//
ProjectionImage::ProjectionImage(Projector::CameraParameters cameraParams, Projector::DetectorParameters imageParams) :
	m_cameraParams{ cameraParams },
	m_imageParams{ imageParams},
	m_annData(nullptr),
	m_edgeFilterWidth(1.6)
{
	size_t arraySize = (size_t) m_imageParams.size.width() * m_imageParams.size.height();
	try {
		m_annData = new AnnotationData[arraySize];
	}
	catch (const std::bad_alloc& e) {
		std::cout << "Distance volume allocation failed: " << e.what() << '\n';
	}
	clear(0);
}
ProjectionImage::ProjectionImage(const char* filename) :
	m_cameraParams{ 0, QMatrix4x4() },
	m_imageParams{ QSize(), QVector2D() },
	m_annData(nullptr),
	m_edgeFilterWidth(1.6)
{
	importAnnotation(filename);
}
ProjectionImage::~ProjectionImage()
{
	delete[] m_annData;
	m_annData = nullptr;
}
void ProjectionImage::clear(float maxDist)
{
	size_t arraySize = (size_t)m_imageParams.size.width() * m_imageParams.size.height();
	AnnotationData* pData = m_annData;
	for (int i = 0; i < arraySize; i++) {
		pData->dist = maxDist;
		pData->dir = QVector2D(1, 0); // Arbitrary vector of unit length
		pData->radius = 0;
		pData->time = 0;
		pData->velocity = QVector2D(0, 0);
		pData++;
	}
}

ProjectionImage::AnnotationData ProjectionImage::dataAtPixel(QVector2D pixel)
{ 
	// Perform bounds checking. Clamp to image boundary. Assumes image width and height >= 2 pixels.
	int width = m_imageParams.size.width();
	int height = m_imageParams.size.height();
	if (pixel.x() < 0) {
		pixel.setX(0);
	}
	if (pixel.x() > width - 2) {
		pixel.setX(width - 2);
	}
	if (pixel.y() < 0) {
		pixel.setY(0);
	}
	if (pixel.y() > height - 2) {
		pixel.setY(height - 2);
	}

	// Compute bilinear interpolation parameters
	int i0 = (int)pixel.x();
	int j0 = (int)pixel.y();
	float dx = pixel.x() - i0;
	float dy = pixel.y() - j0;
	int i1 = i0 + 1;
	int j1 = j0 + 1;
	float ddx = 1 - dx;;
	float ddy = 1 - dy;

	// Perform binlinear interpolation
	AnnotationData& d00 = m_annData[i0 + j0 * width];
	AnnotationData& d01 = m_annData[i0 + j1 * width];
	AnnotationData& d10 = m_annData[i1 + j0 * width];
	AnnotationData& d11 = m_annData[i1 + j1 * width];
	AnnotationData data;
	QVector2D vecDist = ddx * ddy * d00.dist * d00.dir + ddx * dy * d01.dist * d01.dir + 
		dx * ddy * d10.dist * d10.dir + dx * dy * d11.dist * d11.dir;
	data.dist = vecDist.length();
	data.dir = vecDist.normalized();
	data.radius = ddx * ddy * d00.radius + ddx * dy * d01.radius + dx * ddy * d10.radius + dx * dy * d11.radius;
	data.time = ddx * ddy * d00.time + ddx * dy * d01.time + dx * ddy * d10.time + dx * dy * d11.time;
	data.velocity = ddx * ddy * d00.velocity + ddx * dy * d01.velocity + dx * ddy * d10.velocity + dx * dy * d11.velocity;
	return data;
}
float ProjectionImage::distAtPixel(QVector2D pixel)
{
	// Perform bounds checking. Clamp to image boundary. Assumes image width and height >= 2 pixels.
	int width = m_imageParams.size.width();
	int height = m_imageParams.size.height();
	if (pixel.x() < 0) {
		pixel.setX(0);
	}
	if (pixel.x() > width - 2) {
		pixel.setX(width - 2);
	}
	if (pixel.y() < 0) {
		pixel.setY(0);
	}
	if (pixel.y() > height - 2) {
		pixel.setY(height - 2);
	}

	// Compute bilinear interpolation parameters
	int i0 = (int)pixel.x();
	int j0 = (int)pixel.y();
	float dx = pixel.x() - i0;
	float dy = pixel.y() - j0;
	int i1 = i0 + 1;
	int j1 = j0 + 1;
	float ddx = 1 - dx;;
	float ddy = 1 - dy;

	// Perform binlinear interpolation
	AnnotationData& d00 = m_annData[i0 + j0 * width];
	AnnotationData& d01 = m_annData[i0 + j1 * width];
	AnnotationData& d10 = m_annData[i1 + j0 * width];
	AnnotationData& d11 = m_annData[i1 + j1 * width];
	AnnotationData data;
	QVector2D vecDist = ddx * ddy * d00.dist * d00.dir + ddx * dy * d01.dist * d01.dir +
		dx * ddy * d10.dist * d10.dir + dx * dy * d11.dist * d11.dir;
	return vecDist.length();
}

void ProjectionImage::exportProjImages(const char* basename)
{
	Stats stats = annotationStats();
	std::string s = std::string(basename);
	exportDistToCenterlineImage((s + "distToCenterline.jpg").c_str(), stats.maxDist);
	exportDistToVesselImage((s + "distToVessel.jpg").c_str(), stats.maxDist);
	exportRadiusImage((s + "radius.jpg").c_str(), stats.maxRadius);
	exportTimeImage((s + "time.jpg").c_str(), stats.maxTime);
	exportSpeedImage((s + "sharedSpeedComponent.jpg").c_str(), stats.maxSpeed);
}
void ProjectionImage::exportAnnotation(const char* filename)
{
	// Create the output stream 
	std::ofstream fstream(filename, std::ios::binary);
	try {
		if (!fstream) {
			throw std::runtime_error("Cannot create the annotation file.");
		}
	}
	catch (const std::exception& e) {
		std::cout << "Exception " << e.what() << std::endl;
		return;
	}

	// Write the annotation version number and parameters
	std::string version = m_version;
	unsigned int stringLength = version.length();
	fstream.write((char*)&stringLength, sizeof(unsigned int));
	fstream.write(version.c_str(), stringLength);
	fstream.write(reinterpret_cast<char*>(&m_cameraParams), sizeof(m_cameraParams));
	fstream.write(reinterpret_cast<char*>(&m_imageParams), sizeof(m_imageParams));

	// Write the annotation data 
	int dataSizeInBytes = m_imageParams.size.width() * m_imageParams.size.height() * sizeof(AnnotationData);
	fstream.write(reinterpret_cast<char*>(m_annData), dataSizeInBytes);
}


// 
// Private
//
ProjectionImage::Stats ProjectionImage::annotationStats()
{
	AnnotationData* pData = m_annData;
	float maxDist = 0;
	float maxRadius = 0;
	float maxTime = 0;
	float maxSpeed = 0;
	size_t arraySize = (size_t)m_imageParams.size.width() * m_imageParams.size.height();
	for (int i = 0; i < arraySize; i++) {
		if (pData->dist > maxDist) maxDist = pData->dist;
		if (pData->radius > maxRadius) maxRadius = pData->radius;
		if (pData->time > maxTime) maxTime = pData->time;
		float speed = pData->velocity.length();
		if (speed > maxSpeed) {
			maxSpeed = speed;
		}
		pData++;
	}
	return Stats({ maxDist, maxRadius, maxTime, maxSpeed });
}

void ProjectionImage::exportDistToCenterlineImage(const char* filename, float maxDist)
{
	QImage image(m_imageParams.size, QImage::Format_Grayscale8);
	AnnotationData* pData = m_annData;
	for (int j = 0; j < m_imageParams.size.height(); ++j) {
		unsigned char* pImage = image.scanLine(j);
		for (int i = 0; i < m_imageParams.size.width(); ++i) {
			float scaled = std::max(0.0f, std::min(1.0f, (maxDist - pData->dist) / maxDist));
			*pImage++ = (unsigned char)(255.0 * scaled);
			pData++;
		}
	}
	exportImage(filename, image);
}
void ProjectionImage::exportDistToVesselImage(const char* filename, float maxDist)
{
	QImage image(m_imageParams.size, QImage::Format_Grayscale8);
	AnnotationData* pData = m_annData;
	for (int j = 0; j < m_imageParams.size.height(); ++j) {
		unsigned char* pImage = image.scanLine(j);
		for (int i = 0; i < m_imageParams.size.width(); ++i) {
			float distToEdge = pData->radius - pData->dist;
			float scaledOffsetDistToEdge = std::max(0.0f, std::min(1.0f, 0.5f * distToEdge / maxDist + 0.5f));
			*pImage++ = (unsigned char)(255.0 * scaledOffsetDistToEdge);
			pData++;
		}
	}
	exportImage(filename, image);
}
void ProjectionImage::exportRadiusImage(const char* filename, float maxRadius)
{
	QImage image(m_imageParams.size, QImage::Format_Grayscale8);
	AnnotationData* pData = m_annData;
	for (int j = 0; j < m_imageParams.size.height(); ++j) {
		unsigned char* pImage = image.scanLine(j);
		for (int i = 0; i < m_imageParams.size.width(); ++i) {
			float offsetDistToEdge = pData->radius - pData->dist + m_edgeFilterWidth / 2.0;
			float edgeFilter = std::max(0.0f, std::min(1.0f, offsetDistToEdge / m_edgeFilterWidth));
			float scaledRadius = std::max(0.0f, std::min(1.0f, pData->radius / maxRadius));
			*pImage++ = (unsigned char)(255.0f * (scaledRadius * edgeFilter));
			pData++;
		}
	}
	exportImage(filename, image);
}
void ProjectionImage::exportTimeImage(const char* filename, float maxTime)
{
	QImage image(m_imageParams.size, QImage::Format_Grayscale8);
	AnnotationData* pData = m_annData;
	for (int j = 0; j < m_imageParams.size.height(); ++j) {
		unsigned char* pImage = image.scanLine(j);
		for (int i = 0; i < m_imageParams.size.width(); ++i) {
			float offsetDistToEdge = pData->radius - pData->dist + m_edgeFilterWidth / 2.0;
			float edgeFilter = std::max(0.0f, std::min(1.0f, offsetDistToEdge / m_edgeFilterWidth));
			float scaledTime = std::max(0.0f, std::min(1.0f, pData->time / maxTime));
			*pImage++ = (unsigned char)(255.0f * (scaledTime * edgeFilter));
			pData++;
		}
	}
	exportImage(filename, image);
}
void ProjectionImage::exportSpeedImage(const char* filename, float maxSpeed)
{
	QImage image(m_imageParams.size, QImage::Format_Grayscale8);
	AnnotationData* pData = m_annData;
	for (int j = 0; j < m_imageParams.size.height(); ++j) {
		unsigned char* pImage = image.scanLine(j);
		for (int i = 0; i < m_imageParams.size.width(); ++i) {
			float offsetDistToEdge = pData->radius - pData->dist + m_edgeFilterWidth / 2.0;
			float edgeFilter = std::max(0.0f, std::min(1.0f, offsetDistToEdge / m_edgeFilterWidth));
			float signedDyByDt = pData->velocity[1];   // TODO **************************** assumes y direction is common to images and velocity is normalized to 1.0
			float scaledOffsetDyDt = std::max(0.0f, std::min(1.0f, 0.5f * signedDyByDt / maxSpeed + 0.5f));
			*pImage++ = (unsigned char)(255.0f * scaledOffsetDyDt * edgeFilter);
			pData++;
		}
	}
	exportImage(filename, image);
}

void ProjectionImage::exportImage(const char* filename, QImage& image) {
	try {
		if (!image.save(filename, nullptr, 100)) {
			throw std::runtime_error("Can't save image to export file.");
		}
	}
	catch (const std::exception& e) {
		std::cout << "Exception " << e.what() << std::endl;
	}
}
void ProjectionImage::importAnnotation(const char* filename)
{
	// Create the input stream 
	std::ifstream fstream(filename, std::ios::binary);
	try {
		if (!fstream) {
			throw std::runtime_error("Cannot open the annotation file.");
		}
	}
	catch (const std::exception& e) {
		std::cout << "Exception " << e.what() << std::endl;
		return;
	}

	// Read the annotation version number and parameters
	unsigned int stringLength;
	fstream.read((char*)&stringLength, sizeof(unsigned int));
	std::string version;
	version.resize(stringLength);
	fstream.read((char*)version.c_str(), stringLength);
	if (version != m_version) {
		std::cout << "Error. Bad input file or file version." << std::endl;
		return;
	}
	fstream.read(reinterpret_cast<char*>(&m_cameraParams), sizeof(m_cameraParams));
	fstream.read(reinterpret_cast<char*>(&m_imageParams), sizeof(m_imageParams));

	// Allocate memory and read the annotation data 
	if (m_annData) {
		delete[] m_annData;
		m_annData = nullptr;
	}
	size_t arraySize = (size_t)m_imageParams.size.width() * m_imageParams.size.height();
	m_annData = new AnnotationData[arraySize];
	int dataSizeInBytes = arraySize * sizeof(AnnotationData);
	fstream.read(reinterpret_cast<char*>(m_annData), dataSizeInBytes);

	Stats stats = annotationStats();

}
