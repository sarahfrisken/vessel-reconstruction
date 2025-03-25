//
// XmlUtils.h
// Utilities for reading and writing data required by the application to xml files.
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

#include <string>
#include <vector>

#include "Projector.h"
#include "Reconstructor.h"

// TODO, separate into 3 sets of xml utilities: generationXML, reconstructionXML, validationXML,
// each with writeTemplate, read, data access (reconstruction also has write)
class XmlUtils
{
public:
	XmlUtils();
	~XmlUtils();

	static int writeGenerateFromBravaTemplate(std::string filename);
	static int writeReconstructionTemplate(std::string filename);
	static int writeValidationTemplate(std::string filename);

	static int XmlUtils::writeReconstruction(std::string filename, std::vector<std::string> imageFilenames,
		ReconstructionParameters reconParams, ReconstructionThresholds reconThresholds);

	int readSyntheticGenerationXml(std::string filename);
	int readReconstructionXml(std::string filename);
	int readValidationXml(std::string filename);

	// Access state read from a synthetic generation xml file
	std::string modelFilename();
	std::vector<int> circulations();
	int maxDim();
	float maxDistInImage();
	bool doExportVolumes();
	bool doExportAllImages();
	std::vector<Projector::CameraParameters> genCameraParams();
	std::string genExportDir();

	// Access state read from a reconstruction xml file
	std::vector<std::string> imageFilenames();
	ReconstructionParameters reconParams();
	ReconstructionThresholds reconThresholds();
	std::string reconExportDir();

	std::vector<ReconstructionThresholds> validationThresholds();
	std::string validateExportFilename();

private:
	// State read for synthetic generation
	bool m_syntheticGenerationRead;
	bool m_doExportVols;
	bool m_doExportAllImages;
	std::string m_modelFilepath;
	int m_maxDim;
	float m_maxDistInImage;
	std::vector<int> m_circulations;
	std::string m_genExportPath;
	struct CameraData {
		float focalLength = 0;
		QVector3D cameraDir;
		QVector3D cameraUp;
	};
	std::vector<CameraData> m_cameraData;
	QMatrix4x4 worldToCameraMatrix(const QVector3D& cameraPos,
		const QVector3D& targetPos, const QVector3D& desiredUp);

	// State read for reconstruction
	bool m_reconstructionRead;
	std::vector<std::string> m_reconImageFilenames;
	ReconstructionParameters m_reconParams;
	ReconstructionThresholds m_reconThresholds;
	std::string m_reconExportPath;

	// State read only for validation
	struct ThresholdRange {
		float min;
		float max;
		int numSamples;
	};

	bool m_validationRead;
	std::vector<ReconstructionThresholds> m_validationThresholds;
	std::string m_valExportPath;
	std::vector<ReconstructionThresholds> thresholdsFromRanges(ThresholdRange dist,
		ThresholdRange radius, ThresholdRange time, ThresholdRange speed);

	void removeDoubleBackslash(std::string& string);
};