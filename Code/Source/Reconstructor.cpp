//
// Reconstructor.cpp
// Implementation of Reconstructor.
//

#include "Reconstructor.h"

#include <string>
#include <iostream>
#include <fstream>

// 
// Public
//
Reconstructor::Reconstructor(std::vector<std::shared_ptr<ProjectionImage>> projImages, ReconstructionParameters reconParams) :
	m_projImages(projImages),
	m_reconParams(reconParams)
{
	if (m_projImages.size() == 0) return;
	testProjection();
}
Reconstructor::~Reconstructor() 
{
}

void Reconstructor::reconstructAndExport(const char* filename)
{
	// Use default reconstruction thresholds
	ReconstructionThresholds rt;
	reconstructAndExport(filename, rt);
}
void Reconstructor::reconstructAndExport(const char* filename, ReconstructionThresholds thresholds)
{
	switch (m_reconParams.type()) {
	case (ReconstructionParameters::ReconstructionType::MaxDist):
	{
		reconstructAndExportMaxDist(filename);
		break;
	}
	case (ReconstructionParameters::ReconstructionType::AvgDist):
	{
		reconstructAndExportAvgDist(filename);
		break;
	}
	case (ReconstructionParameters::ReconstructionType::ConstrainedDist):
	{
		reconstructAndExportConstrainedDist(filename, thresholds);
		break;
	}
	}
}
std::vector<ValidationMetrics> Reconstructor::reconstructAndValidate(std::shared_ptr<Volume> distField,
	std::vector<ReconstructionThresholds> thresholds)
{
	int numProjections = m_projImages.size();
	if (numProjections == 2) {
		return reconstructAndValidateConstrainedDist(distField, thresholds);
	}
	else {
		ValidationMetrics m = reconstructAndValidateMaxDist(distField, thresholds[0].maxDist());
		std::vector<ValidationMetrics> metrics;
		metrics.push_back(m);
		return metrics;
	}
}


//
// Private
//
QVector3D voxelToWorld(std::array<int, 3> size, float voxelSize, std::array<int, 3> voxel)
{
	// Assumes the reconstruction volume is centered at the world origin and aligned with 
	// the principle axes
	QVector3D pWorld;
	pWorld[0] = (voxel[0] - 0.5f * size[0]) * voxelSize;
	pWorld[1] = (voxel[1] - 0.5f * size[1]) * voxelSize;
	pWorld[2] = (voxel[2] - 0.5f * size[2]) * voxelSize;
	return pWorld;
}
QVector3D voxelToWorld(std::array<int, 3> size, float voxelSize, QVector3D volumeOrigin, 
	std::array<int, 3> voxel)
{
	// Assumes the volume is aligned with the world x, y, z axes and its left, front, 
	// bottom corner is located at volumeOrigin
	// Not used -- instead assume volume is centered at world origin. volumeOrigin is 
	// used with Brava data to encode offset build into vessel endpoints in vesselness
	// volume.
	QVector3D pWorld;
	pWorld *= voxelSize;
	pWorld += volumeOrigin;
	return pWorld;
}
void Reconstructor::reconstructAndExportAvgDist(const char* filename)
{
	Volume vesselnessVol(m_reconParams.volumeSize(), m_reconParams.voxelSize(), m_reconParams.volumeOrigin(), 0);
	
	int numProjections = m_projImages.size();
	if (numProjections < 1) return;

	std::vector<Projector> projectors;
	for (int i = 0; i < numProjections; i++) {
		projectors.push_back(Projector(m_projImages[i]->cameraParameters(), m_projImages[i]->imageParameters()));
	}

	// Perform the reconstruction, storing the average distance in each voxel
	std::array<int, 3> size = m_reconParams.volumeSize();
	float voxelSize = m_reconParams.voxelSize();
	float maxDist = m_reconParams.maxDist();
	float averagingScale = 1.0f / ((float)numProjections);
	for (int k = 0; k < size[2]; k++) {
		for (int j = 0; j < size[1]; j++) {
			for (int i = 0; i < size[0]; i++) {
				std::array<int, 3> voxel{ i, j, k };
				QVector3D worldPoint = voxelToWorld(size, voxelSize, voxel);
				for (int idx = 0; idx < numProjections; idx++) {
					float pixelToWorldLengthScale = projectors[idx].pixelLengthToWorldScale(worldPoint);
					QVector2D pixel = projectors[idx].worldToPixel(worldPoint);
					float dist = m_projImages[idx]->distAtPixel(pixel);
					float reconDist = dist * pixelToWorldLengthScale;
					float vesselness = vesselness = 1 - std::max(0.0f, std::min(1.0f, reconDist / maxDist));
					float stored;
					if (vesselnessVol.dataAtVoxel(voxel, stored)) {
						vesselnessVol.setDataAtVoxel(voxel, stored + vesselness * averagingScale);
					}
				}
			}
		}
	}
	vesselnessVol.exportAsScalarNRRD(filename, 0, 1, false);
}
void Reconstructor::reconstructAndExportMaxDist(const char* filename)
{
	Volume vesselnessVol(m_reconParams.volumeSize(), m_reconParams.voxelSize(), m_reconParams.volumeOrigin(), 1.0);

	int numProjections = m_projImages.size();
	if (numProjections < 1) return;

	std::vector<Projector> projectors;
	for (int i = 0; i < numProjections; i++) {
		projectors.push_back(Projector(m_projImages[i]->cameraParameters(), m_projImages[i]->imageParameters()));
	}
	
	// Perform the reconstruction, storing the average distance in each voxel
	std::array<int, 3> size = m_reconParams.volumeSize();
	float voxelSize = m_reconParams.voxelSize();
	float maxDist = m_reconParams.maxDist();
	for (int k = 0; k < size[2]; k++) {
		for (int j = 0; j < size[1]; j++) {
			for (int i = 0; i < size[0]; i++) {
				std::array<int, 3> voxel{ i, j, k };
				QVector3D worldPoint = voxelToWorld(size, voxelSize, voxel);

				for (int idx = 0; idx < numProjections; idx++) {
					float pixelToWorldLengthScale = projectors[idx].pixelLengthToWorldScale(worldPoint);
					QVector2D pixel = projectors[idx].worldToPixel(worldPoint);
					float dist = m_projImages[idx]->distAtPixel(pixel);
					float reconDist = dist * pixelToWorldLengthScale;
					float vesselness = 1 - std::max(0.0f, std::min(1.0f, reconDist / maxDist));
					float stored;
					if (vesselnessVol.dataAtVoxel(voxel, stored)) {
						if (vesselness < stored) {
							vesselnessVol.setDataAtVoxel(voxel, vesselness);
						}
					}
				}
			}
		}
	}
	vesselnessVol.exportAsScalarNRRD(filename, 0, 1, false);
}
float lerp(float a, float b, float t) {
	// Linear interpolation with clamped output
	if (t < 0.0f) return a;
	if (t > 1.0f) return b;
	return (1.0f - t) * a + t * b;
}
void Reconstructor::reconstructAndExportConstrainedDist(const char* filename, ReconstructionThresholds thresholds)
{
	Volume vesselnessVol(m_reconParams.volumeSize(), m_reconParams.voxelSize(), m_reconParams.volumeOrigin(), 1.0);

	// Handles two projectors. To handle more projectors, use variances rather than differences
	int numProjections = m_projImages.size();
	if (numProjections != 2) return;
	std::vector<Projector> projectors;
	for (int i = 0; i < numProjections; i++) {
		projectors.push_back(Projector(m_projImages[i]->cameraParameters(), m_projImages[i]->imageParameters()));
	}

	// Perform the reconstruction
	std::array<int, 3> size = m_reconParams.volumeSize();
	float voxelSize = m_reconParams.voxelSize();
	float maxDist = m_reconParams.maxDist();
	for (int k = 0; k < size[2]; k++) {
		for (int j = 0; j < size[1]; j++) {
			for (int i = 0; i < size[0]; i++) {
				std::array<int, 3> voxel{ i, j, k };
				QVector3D worldPoint = voxelToWorld(size, voxelSize, voxel);

				ProjectionImage::AnnotationData data0 = m_projImages[0]->dataAtPixel(projectors[0].worldToPixel(worldPoint));
				ProjectionImage::AnnotationData data1 = m_projImages[1]->dataAtPixel(projectors[1].worldToPixel(worldPoint));

				// Back project data to voxel coordinates
				float pixelToWorldLengthScale0 = projectors[0].pixelLengthToWorldScale(worldPoint);
				float pixelToWorldLengthScale1 = projectors[1].pixelLengthToWorldScale(worldPoint);
				float dist0 = data0.dist * pixelToWorldLengthScale0;
				float dist1 = data1.dist * pixelToWorldLengthScale1;
				float radius0 = data0.radius * pixelToWorldLengthScale0;
				float radius1 = data1.radius * pixelToWorldLengthScale1;
				QVector2D velocity0 = data0.velocity * pixelToWorldLengthScale0;
				QVector2D velocity1 = data1.velocity * pixelToWorldLengthScale1;

				float reconDist = std::max(dist0, dist1);
				float vesselness = 1 - std::max(0.0f, std::min(1.0f, reconDist / maxDist));

				// Confidence weights in [0, 1]. Max 1 --> can take vesselness to 0. 
				double w = 0.5;

				// Compute confidence. Confidence is highest when differences are small.
				float diffRadius = fabs(radius0 - radius1);
				float diffTime = fabs(data0.time - data1.time);
				float diffSharedVelocityComponent = fabs(velocity0.y() - velocity1.y());  // ************* TODO assumes y vector is common
				float cRadius = 1.0 - w * std::min(1.0f, diffRadius / thresholds.maxRadiusDiff());
				float cTime = 1.0 - w * std::min(1.0f, diffTime / thresholds.maxTimeDiff());
				float cSpeed = 1.0 - w * std::min(1.0f, diffSharedVelocityComponent / thresholds.maxSharedVelCompDiff());
				float confidence = cRadius * cTime * cSpeed;
				vesselness *= confidence;
				vesselnessVol.setDataAtVoxel(voxel, vesselness);
			}
		}
	}
	vesselnessVol.exportAsScalarNRRD(filename, 0, 1, false);
}
std::vector<ValidationMetrics> Reconstructor::reconstructAndValidateConstrainedDist(std::shared_ptr<Volume> distField,
	std::vector<ReconstructionThresholds> thresholds)
{
	std::vector<ValidationMetrics> metrics;
	float* pDistFieldData = distField->data();
	if (!pDistFieldData) return metrics;

	// Handles two projectors
	int numProjections = m_projImages.size();
	if (numProjections != 2) return metrics;
	std::vector<Projector> projectors;
	for (int i = 0; i < numProjections; i++) {
		projectors.push_back(Projector(m_projImages[i]->cameraParameters(), m_projImages[i]->imageParameters()));
	}

	// Initialize the validation metrics
	for (int i = 0; i < thresholds.size(); i++) {
		ValidationMetrics m;
		metrics.push_back(m);
	}

	// Perform the validation. Reconstruct the constrained distance at each voxel and compute 
	// validation metrics based on distance field at that voxel.
	std::array<int, 3> size = m_reconParams.volumeSize();
	float voxelSize = m_reconParams.voxelSize();
	float maxDist = m_reconParams.maxDist();
	for (int k = 0; k < size[2]; k++) {
		for (int j = 0; j < size[1]; j++) {
			for (int i = 0; i < size[0]; i++) {
				std::array<int, 3> voxel{ i, j, k };
				QVector3D worldPoint = voxelToWorld(size, voxelSize, voxel);

				ProjectionImage::AnnotationData data0 = m_projImages[0]->dataAtPixel(projectors[0].worldToPixel(worldPoint));
				ProjectionImage::AnnotationData data1 = m_projImages[1]->dataAtPixel(projectors[1].worldToPixel(worldPoint));

				float distAtPoint(0);
				bool inBounds = distField->dataAtWorldPoint(worldPoint, distAtPoint);
				if (!inBounds) {
					// Do not count this voxel. Do not update the validation metrics.
					continue;
				}

				// Back project data to voxel coordinates
				float pixelToWorldLengthScale0 = projectors[0].pixelLengthToWorldScale(worldPoint);
				float pixelToWorldLengthScale1 = projectors[1].pixelLengthToWorldScale(worldPoint);
				float dist0 = data0.dist * pixelToWorldLengthScale0;
				float dist1 = data1.dist * pixelToWorldLengthScale1;
				float radius0 = data0.radius * pixelToWorldLengthScale0;
				float radius1 = data1.radius * pixelToWorldLengthScale1;
				QVector2D velocity0 = data0.velocity * pixelToWorldLengthScale0;
				QVector2D velocity1 = data1.velocity * pixelToWorldLengthScale1;

				float reconDist = std::max(dist0, dist1);
				float vesselness = 1 - std::max(0.0f, std::min(1.0f, reconDist / maxDist));

				// Compute confidence. Confidence is highest when differences are small.
				float diffRadius = fabs(radius0 - radius1);
				float diffTime = fabs(data0.time - data1.time);
				float diffSharedVelocityComponent = fabs(velocity0.y() - velocity1.y());  // ************* TODO assumes y vector is common
		
				// Perform validation for each set of thresholds
				for (int idx = 0; idx < thresholds.size(); idx++) {
					bool doHardTest = false;
					bool passAll = true;
					if (doHardTest) {
						if (reconDist >= thresholds[idx].maxDist()) passAll = false;
						else if (diffRadius >= thresholds[idx].maxRadiusDiff()) passAll = false;
						else if (diffTime >= thresholds[idx].maxTimeDiff()) passAll = false;
						else if (diffSharedVelocityComponent >= thresholds[idx].maxSharedVelCompDiff()) passAll = false;
					}
					else {
						float cRadius = 1.0 - std::min(1.0f, diffRadius / thresholds[idx].maxRadiusDiff());
						float cTime = 1.0 - std::min(1.0f, diffTime / thresholds[idx].maxTimeDiff());
						float cSpeed = 1.0 - std::min(1.0f, diffSharedVelocityComponent / thresholds[idx].maxSharedVelCompDiff());
						float confidence = cRadius * cTime * cSpeed;
						float modifiedVesselness = vesselness * confidence;
						passAll = ((1 - modifiedVesselness) * maxDist < thresholds[idx].maxDist()) ? true : false;
					}

					if (distAtPoint < thresholds[idx].maxDist()) {
						// Point is in vessel
						if (passAll) metrics[idx].incrementTruePos();
						else metrics[idx].incrementFalseNeg();
					}
					else {
						if (passAll) metrics[idx].incrementFalsePos();
						else metrics[idx].incrementTrueNeg();
					}
				}
			}
		}
	}

	return metrics;
}
ValidationMetrics Reconstructor::reconstructAndValidateMaxDist(std::shared_ptr<Volume> distField, float maxDistThreshold)
{
	ValidationMetrics metrics;
	float* pDistFieldData = distField->data();
	if (!pDistFieldData) return metrics;

	// Handles two projectors
	int numProjections = m_projImages.size();
	if (numProjections < 1) return metrics;
	std::vector<Projector> projectors;
	for (int i = 0; i < numProjections; i++) {
		projectors.push_back(Projector(m_projImages[i]->cameraParameters(), m_projImages[i]->imageParameters()));
	}

	// Perform the validation. Reconstruct the distance at each voxel and compute validation metrics based on 
	// distance field at that voxel.
	std::array<int, 3> size = m_reconParams.volumeSize();
	float voxelSize = m_reconParams.voxelSize();
	float maxDist = m_reconParams.maxDist();
	for (int k = 0; k < size[2]; k++) {
		for (int j = 0; j < size[1]; j++) {
			for (int i = 0; i < size[0]; i++) {
				std::array<int, 3> voxel{ i, j, k };
				QVector3D worldPoint = voxelToWorld(size, voxelSize, voxel);

				float distAtPoint(0);
				bool inBounds = distField->dataAtWorldPoint(worldPoint, distAtPoint);
				if (!inBounds) {
					// Do not count this voxel. Do not update the validation metrics.
					continue;
				}

				float maxDist = 0;
				for (int idx = 0; idx < numProjections; idx++) {
					float pixelToWorldLengthScale = projectors[idx].pixelLengthToWorldScale(worldPoint);
					float dist = m_projImages[idx]->distAtPixel(projectors[idx].worldToPixel(worldPoint));
					float reconDist = dist * pixelToWorldLengthScale;
					if (reconDist > maxDist) maxDist = reconDist;
				}

				if (distAtPoint < maxDistThreshold) {
					// Point is in vessel
					if (maxDist < maxDistThreshold) metrics.incrementTruePos();
					else metrics.incrementFalseNeg();
				}
				else {
					if (maxDist < maxDistThreshold) metrics.incrementFalsePos();
					else metrics.incrementTrueNeg();
				}
			}
		}
	}
	return metrics;
}

void Reconstructor::testProjection()
{
	std::vector<Projector> projectors;
	int numProjections = m_projImages.size();
	for (int i = 0; i < numProjections; i++) {
		projectors.push_back(Projector(m_projImages[i]->cameraParameters(), m_projImages[i]->imageParameters()));
	}

	std::array<int, 3> size = m_reconParams.volumeSize();
	float voxelSize = m_reconParams.voxelSize();
	std::array<int, 3> voxel{ size[0] / 2, size[1] / 2, size[2] / 2 };
	QVector3D volumeCenterInWorldCoords = voxelToWorld(size, voxelSize, voxel);

	for (int idx = 0; idx < numProjections; idx++) {
		QVector4D imagePoint = projectors[idx].worldToImage(volumeCenterInWorldCoords);
		QVector2D pixel = projectors[idx].worldToPixel(volumeCenterInWorldCoords);

		QVector4D cameraCenterInWorldCoords = projectors[idx].cameraToWorld(QVector3D(0, 0, 750));
		imagePoint = projectors[idx].worldToImage(QVector3D(cameraCenterInWorldCoords));
		pixel = projectors[idx].worldToPixel(QVector3D(cameraCenterInWorldCoords));

		QVector4D cameraPointInWorldCoords = projectors[idx].cameraToWorld(QVector3D(10, 10, 750));
		imagePoint = projectors[idx].worldToImage(QVector3D(cameraPointInWorldCoords));
		pixel = projectors[idx].worldToPixel(QVector3D(cameraPointInWorldCoords));
		pixel += {0, 0};
	}
}