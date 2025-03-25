//
// VesselModel.cpp
// Implementation of VesselModel.
//

#include "VesselModel.h"

#include <string>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <QVector2D>

// 
// Public
//
VesselModel::VesselModel(std::vector<VesselModel::Line> lines) :
	m_lines(lines),
	m_bbox{ {0, 0, 0}, {0, 0, 0} }
{
	centerAndComputeBBox();
}
VesselModel::~VesselModel()
{
}

VesselModel::BBox VesselModel::boundingBox()
{
	return m_bbox;
}


void VesselModel::exportDistToCenterlineVolume(const char* filename, float voxelSize, 
	float maxDist)
{
	QVector3D offset;
	std::array<int, 3> size = getEnclosingVolumeSize(m_bbox, voxelSize, maxDist, offset);
	Volume vol(size, voxelSize, offset, maxDist);

	// Scan model lines into the volume, storing minimum distance to the closest vessel 
	// centerline
	for (std::vector<VesselModel::Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		addCenterline(vol, *it, maxDist);
	}

	// Export the distance volume 
	bool invert = true;
	vol.exportAsScalarNRRD(filename, 0, maxDist, invert);
}
void VesselModel::exportSignedDistToVesselVolume(const char* filename, float voxelSize, 
	float maxDist)
{
	QVector3D offset;
	std::array<int, 3> size = getEnclosingVolumeSize(m_bbox, voxelSize, maxDist, offset);
	Volume vol(size, voxelSize, offset, maxDist);

	// Scan model lines into the volume, storing minimum distance to the closest vessel 
	// edge. Distances are positive outside the vessel and negative inside.
	for (std::vector<VesselModel::Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		addVessel(vol, *it, maxDist);
	}

	// Export the distance volume 
	bool invert = true;
	vol.exportAsScalarNRRD(filename, -maxDist, maxDist, invert);
}

std::shared_ptr<Volume>  VesselModel::distToCenterlineVolume(float voxelSize, float maxDist)
{
	// Get volume size and offset
	QVector3D offset;
	std::array<int, 3> size = getEnclosingVolumeSize(m_bbox, voxelSize, maxDist, offset);
	std::shared_ptr<Volume> pVol(new Volume(size, voxelSize, offset, maxDist));

	// Scan model lines into the volume, storing minimum distance to the closest vessel 
	// centerline
	for (std::vector<VesselModel::Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		addCenterline(*pVol, *it, maxDist);
	}
	return pVol;
}

Projector::DetectorParameters VesselModel::getOptimalImageParams(
	Projector::CameraParameters cameraParams, float maxDim, float maxDistInImage)
{
	// Project each corner of the model bounding box onto the image
	// detParams are arbitrary as they are not needed by worldToImage()
	Projector::DetectorParameters detParams = {
		QSize(0, 0),
		QVector2D(1.0, 1.0)
	};
	Projector p(cameraParams, detParams);

	QVector3D corner[8];
	corner[0] = QVector3D(m_bbox.min[0], m_bbox.min[1], m_bbox.min[2]);
	corner[1] = QVector3D(m_bbox.max[0], m_bbox.min[1], m_bbox.min[2]);
	corner[2] = QVector3D(m_bbox.max[0], m_bbox.max[1], m_bbox.min[2]);
	corner[3] = QVector3D(m_bbox.min[0], m_bbox.max[1], m_bbox.min[2]);
	corner[4] = QVector3D(m_bbox.min[0], m_bbox.min[1], m_bbox.max[2]);
	corner[5] = QVector3D(m_bbox.max[0], m_bbox.min[1], m_bbox.max[2]);
	corner[6] = QVector3D(m_bbox.max[0], m_bbox.max[1], m_bbox.max[2]);
	corner[7] = QVector3D(m_bbox.min[0], m_bbox.max[1], m_bbox.max[2]);

	QVector2D min(p.worldToImage(corner[0]));
	QVector2D max = min;
	for (int i = 1; i < 8; i++) {
		QVector2D proj = QVector2D(p.worldToImage(corner[i]));
		min[0] = std::min(min[0], proj[0]);
		min[1] = std::min(min[1], proj[1]);
		max[0] = std::max(max[0], proj[0]);
		max[1] = std::max(max[1], proj[1]);
	}
	
	float xSize = max[0] - min[0] + 2 * maxDistInImage;
	float ySize = max[1] - min[1] + 2 * maxDistInImage;
	if (xSize >= ySize) {
		float pixelSize = xSize / maxDim;
		detParams.size = QSize(maxDim, (int) (ySize / pixelSize + 0.5));
		detParams.pixelSize = QVector2D(pixelSize, pixelSize);
	}
	else {
		float pixelSize = ySize / maxDim;
		detParams.size = QSize((int)(xSize / pixelSize + 0.5), maxDim);
		detParams.pixelSize = QVector2D(pixelSize, pixelSize);
	}
	return detParams;
}
void VesselModel::renderProjectionImage(std::shared_ptr<ProjectionImage> image, float maxDistInImage)
{
	// Determine the image parameters, create and initialize the image
	Projector p(image->cameraParameters(), image->imageParameters());

	float maxRadius = maxProjectedRadius(p);	// Distance in mm in the projected image
	float maxOffset = maxDistInImage + maxRadius;		// Offset is in mm in the projected image. It is big enough to enclose the thick line/vessel + the distance field to maxDist
	image->clear(maxOffset);
	
	// Project each line of the model onto the image
	for (std::vector<Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		ProjectedVertex start = projectToPixelCoords(p, it->start);
		ProjectedVertex end = projectToPixelCoords(p, it->end);
		QVector2D lineVec = end.pos - start.pos;
		ProjectedLine projLine = { lineVec.length(), lineVec.normalized(), start, end };
		addLineToProjectedImage(image, projLine, maxOffset);
	}
}

// 
// Private
//
VesselModel::ProjectedVertex VesselModel::projectToPixelCoords(Projector projector,
	VesselModel::Vertex vertex)
{
	QVector2D projPos = projector.worldToPixel(vertex.pos);
	float projRadius = vertex.radius * projector.worldLengthToPixelScale(vertex.pos);
	QVector2D projVelocity = projector.worldVectorToPixel(vertex.pos, vertex.velocity);
	ProjectedVertex projVertex = { projPos, projRadius, vertex.time, projVelocity };
	return projVertex;
}
float VesselModel::maxProjectedRadius(Projector projector)
{
	float maxProjectedRadius = 0;
	for (std::vector<Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		ProjectedVertex vStart = projectToPixelCoords(projector, it->start);
		ProjectedVertex vEnd = projectToPixelCoords(projector, it->end);
		maxProjectedRadius = std::max(maxProjectedRadius, vStart.radius);
		maxProjectedRadius = std::max(maxProjectedRadius, vEnd.radius);
	}
	return maxProjectedRadius;
}
void VesselModel::addLineToProjectedImage(std::shared_ptr<ProjectionImage> image,
	VesselModel::ProjectedLine line, float maxOffset)
{
	ProjectionImage::AnnotationData* pData = image->data();

	// Determine line bounding box
	BBox2D bbox = getProjectedLineBBox(image, line, maxOffset);

	// Update annotated image inside the bounding box
	QSize size = image->imageParameters().size;
	for (int j = bbox.min[1]; j < bbox.max[1]; j++) {
		int jIdx = j * size.width();
		for (int i = bbox.min[0]; i < bbox.max[0]; i++) {
			int idx = i + jIdx;
			QVector2D point { (float)i, (float)j };
			ProjectionImage::AnnotationData data = lineDataAtPoint(point, line);
			if (data.dist < pData[idx].dist) {
				pData[idx] = data;
			}
		}
	}
}
VesselModel::BBox2D VesselModel::getProjectedLineBBox(std::shared_ptr<ProjectionImage> image,
	VesselModel::ProjectedLine line, float maxOffset)
{
	float min[2];
	float max[2];
	for (int i = 0; i < 2; i++) {
		min[i] = std::min(line.start.pos[i], line.end.pos[i]) - maxOffset;
		max[i] = std::max(line.start.pos[i], line.end.pos[i]) + maxOffset;
	}

	// Crop to image
	BBox2D bbox;
	QSize size = image->imageParameters().size;
	bbox.min[0] = std::max(0, (int) min[0]);
	bbox.min[1] = std::max(0, (int) min[1]);
	bbox.max[0] = std::min(size.width() - 1, (int)(max[0] + 0.5));
	bbox.max[1] = std::min(size.height() - 1, (int)(max[1] + 0.5));
	return bbox;
}

ProjectionImage::AnnotationData VesselModel::lineDataAtPoint(QVector2D point,
	VesselModel::ProjectedLine line)
{
	// Find the point on the line closest to point
	float t = 0;
	QVector2D startToP = point - line.start.pos;
	float proj = QVector2D::dotProduct(startToP, line.dir);
	if (proj <= 0 || line.length == 0) t = 0;	// Point is closest to start point
	else if (proj >= line.length) t = 1;		// Point is closest to end point
	else t = proj / line.length;				// Point projects onto the line

	QVector2D closestPoint = (1 - t) * line.start.pos + t * line.end.pos;
	QVector2D vecDistFromClosestPoint = point - closestPoint;

	ProjectionImage::AnnotationData data;
	data.dist = vecDistFromClosestPoint.length();
	data.dir = vecDistFromClosestPoint.normalized();
	data.radius = (1 - t) * line.start.radius + t * line.end.radius;
	data.time = (1 - t) * line.start.time + t * line.end.time;
	data.velocity = (1 - t) * line.start.velocity + t * line.end.velocity;
	return data;
}

void VesselModel::centerAndComputeBBox()
{
	if (m_lines.size() < 1) return;
	std::vector<VesselModel::Line>::iterator it = m_lines.begin();
	BBox bbox;
	bbox.min = m_lines[0].start.pos;
	bbox.max = m_lines[0].start.pos;
	float maxRadius = it->start.radius;
	for (std::vector<VesselModel::Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		// Consider both endpoints of each line to find bounding box
		if (it->start.pos[0] < bbox.min[0]) bbox.min[0] = it->start.pos[0];
		if (it->start.pos[1] < bbox.min[1]) bbox.min[1] = it->start.pos[1];
		if (it->start.pos[2] < bbox.min[2]) bbox.min[2] = it->start.pos[2];
		if (it->end.pos[0] < bbox.min[0]) bbox.min[0] = it->end.pos[0];
		if (it->end.pos[1] < bbox.min[1]) bbox.min[1] = it->end.pos[1];
		if (it->end.pos[2] < bbox.min[2]) bbox.min[2] = it->end.pos[2];

		if (it->start.pos[0] > bbox.max[0]) bbox.max[0] = it->start.pos[0];
		if (it->start.pos[1] > bbox.max[1]) bbox.max[1] = it->start.pos[1];
		if (it->start.pos[2] > bbox.max[2]) bbox.max[2] = it->start.pos[2];
		if (it->end.pos[0] > bbox.max[0]) bbox.max[0] = it->end.pos[0];
		if (it->end.pos[1] > bbox.max[1]) bbox.max[1] = it->end.pos[1];
		if (it->end.pos[2] > bbox.max[2]) bbox.max[2] = it->end.pos[2];

		maxRadius = std::max(maxRadius, it->start.radius);
		maxRadius = std::max(maxRadius, it->end.radius);
	}

	// Center the model
	QVector3D center = 0.5 * (bbox.max + bbox.min);
	for (std::vector<VesselModel::Line>::iterator it = m_lines.begin(); it != m_lines.end(); it++) {
		// Consider both endpoints of each line to find bounding box
		it->start.pos -= center;
		it->end.pos -= center;
	}

	// Incorporate the centering and the vessel radius
	QVector3D radiusOffset(maxRadius, maxRadius, maxRadius);
	m_bbox.min = bbox.min - center - radiusOffset;
	m_bbox.max = bbox.max - center + radiusOffset;
}

std::array<int, 3> VesselModel::getEnclosingVolumeSize(VesselModel::BBox bbox, float voxelSize, float maxDist, QVector3D& offset)
{
	QVector3D paddedMin = { bbox.min[0] - maxDist, bbox.min[1] - maxDist, bbox.min[2] - maxDist };
	QVector3D paddedMax = { bbox.max[0] + maxDist, bbox.max[1] + maxDist, bbox.max[2] + maxDist };
	BBox paddedBBox = { paddedMin, paddedMax };
	QVector3D fSize = paddedBBox.max - paddedBBox.min;
	fSize *= (1.0 / voxelSize);
	std::array<int, 3> size = { int(fSize[0] + 0.5f), int(fSize[1] + 0.5f), int(fSize[2] + 0.5f) };
	offset = paddedBBox.min;
	return size;
}

void VesselModel::addCenterline(Volume& vol, VesselModel::Line& line, float maxDist)
{
	// Convert line endpoint to voxel coordinates
	QVector3D startPos = line.start.pos - vol.origin();
	startPos *= (1.0 / vol.voxelSize());
	QVector3D endPos = line.end.pos - vol.origin();
	endPos *= (1.0 / vol.voxelSize());

	// Determine bounding box of the line and its distance field in voxels. Crop
	// the bounding box to be inside the volume.
	float maxDistInVoxels = maxDist / vol.voxelSize();
	std::array<int, 3> size = vol.size();
	int min[3] = { (int)size[0], (int)size[1], (int)size[2] };
	int max[3] = { 0, 0, 0 };
	for (int i = 0; i < 3; i++) {
		min[i] = std::min(min[i], (int)(std::min(startPos[i], endPos[i]) - maxDistInVoxels));
		min[i] = std::max(0, min[i]);
		max[i] = std::max(max[i], (int)(std::max(startPos[i], endPos[i]) + maxDistInVoxels));
		max[i] = std::min((int)vol.size()[i] - 1, max[i]);
	}

	// Update the distance field around the line
	float* data = vol.data();
	for (int k = min[2]; k < max[2]; k++) {
		int kIdx = k * size[0] * size[1];
		for (int j = min[1]; j < max[1]; j++) {
			int jIdx = j * size[0];
			for (int i = min[0]; i < max[0]; i++) {
				int idx = i + jIdx + kIdx;
				QVector3D p = { (float)i, (float)j, (float)k };
				float distInMM = distPToLine(startPos, endPos, p) * vol.voxelSize();
				if (distInMM < data[idx]) {
					data[idx] = distInMM;
				}
			}
		}
	}
}

void VesselModel::addVessel(Volume& vol, VesselModel::Line& line, float maxDist)
{
	// Convert line endpoint and radius to voxel coordinates
	QVector3D startPos = line.start.pos - vol.origin();
	startPos *= (1.0 / vol.voxelSize());
	QVector3D endPos = line.end.pos - vol.origin();
	endPos *= (1.0 / vol.voxelSize());
	float startRadius = line.start.radius / vol.voxelSize();
	float endRadius = line.end.radius / vol.voxelSize();
	float maxRadius = std::max(startRadius, endRadius);

	// Determine bounding box of the line and its distance field in voxels. Crop
	// the bounding box to be inside the volume.
	float maxDistInVoxels = maxDist / vol.voxelSize();
	std::array<int, 3> size = vol.size();
	int min[3] = { (int)size[0], (int)size[1], (int)size[2] };
	int max[3] = { 0, 0, 0 };
	for (int i = 0; i < 3; i++) {
		min[i] = std::min(min[i], (int)(std::min(startPos[i], endPos[i]) - maxDistInVoxels - maxRadius));
		min[i] = std::max(0, min[i]);
		max[i] = std::max(max[i], (int)(std::max(startPos[i], endPos[i]) + maxDistInVoxels + maxRadius));
		max[i] = std::min((int) vol.size()[i] - 1, max[i]);
	}

	// Update the distance field around the line
	float* data = vol.data();
	for (int k = min[2]; k < max[2]; k++) {
		int kIdx = k * size[0] * size[1];
		for (int j = min[1]; j < max[1]; j++) {
			int jIdx = j * size[0];
			for (int i = min[0]; i < max[0]; i++) {
				int idx = i + jIdx + kIdx;
				QVector3D p = { (float)i, (float)j, (float)k };
				float distInMM = distPToVessel(startPos, startRadius, endPos, endRadius, p) * vol.voxelSize();
				if (distInMM < data[idx]) {
					data[idx] = distInMM;
				}
			}
		}
	}
}

float VesselModel::distPToLine(QVector3D startPos, QVector3D endPos, QVector3D p)
{
	QVector3D pToStart = startPos - p;
	QVector3D lineVec = endPos - startPos;
	float len = lineVec.length();
	if (len == 0) {
		return pToStart.length();
	}
	QVector3D lineDir = lineVec.normalized();
	float proj = -QVector3D::dotProduct(pToStart, lineDir);
	if (proj <= 0) {
		return pToStart.length();
	}
	else if (proj >= len) {
		QVector3D pToEnd = endPos - p;
		return pToEnd.length();
	}
	else {
		float lenPToStart = pToStart.length();
		float distToLineSqr = lenPToStart * lenPToStart - proj * proj;
		if (distToLineSqr <= 0) {
			return 0;
		}
		else {
			return sqrt(distToLineSqr);
		}
	}
}
float VesselModel::distPToVessel(QVector3D startPos, float startRadius, QVector3D endPos, 
	float endRadius, QVector3D p)
{
	QVector3D pToStart = startPos - p;
	QVector3D lineVec = endPos - startPos;
	float len = lineVec.length();
	if (len == 0) {
		return pToStart.length() - startRadius;
	}
	QVector3D lineDir = lineVec.normalized();
	float proj = -QVector3D::dotProduct(pToStart, lineDir);
	if (proj <= 0) {
		return pToStart.length() - startRadius;
	}
	else if (proj >= len) {
		QVector3D pToEnd = endPos - p;
		return pToEnd.length() - endRadius;
	}
	else {
		double t = proj / len;
		float radius = (1.0 - t) * startRadius + t * endRadius;
		float lenPToStart = pToStart.length();
		float distToLineSqr = lenPToStart * lenPToStart - proj * proj;
		if (distToLineSqr <= 0) {
			return 0 - radius;
		}
		else {
			return sqrt(distToLineSqr) - radius;
		}
	}
}