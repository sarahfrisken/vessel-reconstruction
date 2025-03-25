//
// ProjectionImage.h
// Representation and manipulation of annotated projected 2D images of vessels
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

#include "Projector.h"
#include <array>

#include <QImage>
#include <QVector2D>
#include <QVector3D>

class ProjectionImage
{
public:
	ProjectionImage(Projector::CameraParameters cameraParams, Projector::DetectorParameters imageParams);
	ProjectionImage(const char* filename);
	~ProjectionImage();
	void clear(float maxDist);

	struct AnnotationData {
		float dist;			// Distance to nearest vessel centerline point, P
		QVector2D dir;		// Direction to P
		float radius;		// Vessel radius at P
		float time;			// Contrast arrival time at P
		QVector2D velocity;	// Projected contrast velocity at P, (dx/dt, dy/dt)
	};

	ProjectionImage(const ProjectionImage& temp_obj) = delete;
	ProjectionImage& operator=(const ProjectionImage& temp_obj) = delete;

	Projector::CameraParameters cameraParameters() { return m_cameraParams; };
	Projector::DetectorParameters imageParameters() { return m_imageParams; };

	// A pointer to data is provided for efficient data editing and to avoid having 
	// to make copies of the data. External editing must respect the image size to 
	// avoid overflow. If the data is changed, the annotation stats should be updated 
	// before export.
	AnnotationData* data() { return m_annData; };
	AnnotationData dataAtPixel(QVector2D pixel);
	float distAtPixel(QVector2D pixel);

	void exportProjImages(const char* basename);
	void exportAnnotation(const char* basename);

private:
	const std::string m_version = "ANN_IMG_0.1";
	Projector::CameraParameters m_cameraParams;
	Projector::DetectorParameters m_imageParams;
	AnnotationData* m_annData;

	// Annotation stats. Must be explicitly computed using annotationStats()
	struct Stats {
		float maxDist;
		float maxRadius;
		float maxTime;
		float maxSpeed;
	};

	float m_edgeFilterWidth;
	Stats annotationStats();
	void exportDistToCenterlineImage(const char* filename, float maxDist);
	void exportDistToVesselImage(const char* filename, float maxDist);
	void exportRadiusImage(const char* filename, float maxRadius);
	void exportTimeImage(const char* filename, float maxTime);
	void exportSpeedImage(const char* filename, float maxSpeed);
	void exportImage(const char* filename, QImage& image);

	void importAnnotation(const char* filename);
};