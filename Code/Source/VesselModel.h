//
// VesselModel.h
// Represents a 3D vasculature.
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
#include "Projector.h"
#include "ProjectionImage.h"

#include <QVector3D>
#include <vector>

class VesselModel
{
public:
	struct Vertex {
		QVector3D pos;
		float radius;
		float time;
		QVector3D velocity;
	};
	struct Line {
		Vertex start;
		Vertex end;
	};

	VesselModel(std::vector<VesselModel::Line> lines);
	~VesselModel();

	struct BBox {
		QVector3D min;
		QVector3D max;
	};
	BBox boundingBox();
	
	// maxDist is in mm in the model space (i.e., of the exported volume)
	void exportDistToCenterlineVolume(const char* filename, float voxelSize = 1.0, float maxDist = 10.0);
	void exportSignedDistToVesselVolume(const char* filename, float voxelSize = 1.0, float maxDist = 10.0);
	std::shared_ptr<Volume> distToCenterlineVolume(float voxelSize = 1.0, float maxDist = 10.0);

	// Generate a projection image of the vessel model. maxDistInImage is in mm in the projected image.  
	// maxDist in a smaller when backprojected into model space, so scale up accordingly (Scale depends on
	// projection geometry).
	Projector::DetectorParameters getOptimalImageParams(Projector::CameraParameters cameraParams, float maxDim, float maxDistInImage);
	void renderProjectionImage(std::shared_ptr<ProjectionImage> image, float maxDistInImage);

private:
	std::vector<VesselModel::Line> m_lines;

	BBox m_bbox;
	void centerAndComputeBBox();

	// Generate projected images
	struct ProjectedVertex {
		QVector2D pos;
		float radius;
		float time;
		QVector2D velocity;
	};
	struct ProjectedLine {
		float length;
		QVector2D dir;
		ProjectedVertex start;
		ProjectedVertex end;
	};
	struct BBox2D {
		int min[2];
		int max[2];
	};
	ProjectedVertex projectToPixelCoords(Projector projector, Vertex vertex);
	float maxProjectedRadius(Projector projector);
	void addLineToProjectedImage(std::shared_ptr<ProjectionImage> image, ProjectedLine line, float maxOffset);
	BBox2D getProjectedLineBBox(std::shared_ptr<ProjectionImage> image, ProjectedLine line, float maxOffset);
	ProjectionImage::AnnotationData lineDataAtPoint(QVector2D point, ProjectedLine line);

	// Generate distance field volumes
	std::array<int, 3> getEnclosingVolumeSize(BBox bbox, float voxelSize, float maxDist, QVector3D& offset);

	void addCenterline(Volume& vol, VesselModel::Line& line, float maxDist);
	void addVessel(Volume& vol, VesselModel::Line& line, float maxDist);
	float distPToLine(QVector3D startPos, QVector3D endPos, QVector3D p);
	float distPToVessel(QVector3D startPos, float startRadius, QVector3D endPos, float endRadius, QVector3D p);
};