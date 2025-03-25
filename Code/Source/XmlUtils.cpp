//
// XmlUtils.cpp
// Implementation of XmlUtils.
//

#include "XmlUtils.h"

#include "Projector.h"
#include "Reconstructor.h"

#include <string>
#include <vector>
#include "tinyxml2.h"

using namespace tinyxml2;
#ifndef XMLCheckResult
#define XMLCheckResult(a_eResult) if (a_eResult != XML_SUCCESS) { printf("Error: %i\n", a_eResult); return a_eResult; }
#endif

// 
// Public
//
XmlUtils::XmlUtils() :
    m_maxDim(0),
    m_maxDistInImage(10.0f),
    m_syntheticGenerationRead(false),
    m_doExportVols(false),
    m_doExportAllImages(false),
    m_reconstructionRead(false), 
    m_validationRead(false)
{
}
XmlUtils::~XmlUtils()
{
}

// 
// Public
//
int XmlUtils::writeGenerateFromBravaTemplate(std::string filename)
{
    XMLDocument xmlDoc;
    XMLNode* pRoot = xmlDoc.NewElement("Root");
    xmlDoc.InsertEndChild(pRoot);

    XMLElement* pElement = xmlDoc.NewElement("BravaFilepath");
    pElement->SetText("C:/Users/Sarah/Desktop/TestRecon/BG01_ColorCoded.swc");
    pRoot->InsertEndChild(pElement);

    pElement = xmlDoc.NewElement("MaxImageDimension");
    pElement->SetAttribute("maxDim", 256);
    pRoot->InsertEndChild(pElement);

    pElement = xmlDoc.NewElement("MaxDistanceInImage");
    pElement->SetAttribute("maxDist", 3.0);
    pRoot->InsertEndChild(pElement);

    pElement = xmlDoc.NewElement("Circulation");
    pElement->SetAttribute("id", 2);
    pRoot->InsertEndChild(pElement);
    pElement = xmlDoc.NewElement("Circulation");
    pElement->SetAttribute("id", 3);
    pRoot->InsertEndChild(pElement);
    pElement = xmlDoc.NewElement("Circulation");
    pElement->SetAttribute("id", 6);
    pRoot->InsertEndChild(pElement);

    // First projection (AP)
    pElement = xmlDoc.NewElement("CameraData");
    XMLElement* pChild = xmlDoc.NewElement("CameraFocalLength");
    pChild->SetAttribute("f", 1000);    // In mm
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("DirectionToCamera");
    pChild->SetAttribute("x", 0.0);
    pChild->SetAttribute("y", -1.0);    // AP view
    pChild->SetAttribute("z", 0.0);    
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("CameraUpVector");
    pChild->SetAttribute("x", 0.0);
    pChild->SetAttribute("y", 0.0);
    pChild->SetAttribute("z", 1.0);     // IS direction
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Second projection (RL)
    pElement = xmlDoc.NewElement("CameraData");
    pChild = xmlDoc.NewElement("CameraFocalLength");
    pChild->SetAttribute("f", 1000);    // In mm
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("DirectionToCamera");
    pChild->SetAttribute("x", 1.0);
    pChild->SetAttribute("y", 0.0);    // RL view
    pChild->SetAttribute("z", 0.0);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("CameraUpVector");
    pChild->SetAttribute("x", 0.0);
    pChild->SetAttribute("y", 0.0);
    pChild->SetAttribute("z", 1.0);     // IS direction
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Export options
    pElement = xmlDoc.NewElement("ExportOptions");
    pElement->SetAttribute("doExportVolumes", 1);
    pElement->SetAttribute("doExportAllImages", 1);
    pRoot->InsertEndChild(pElement);

    // Exporting directory
    pElement = xmlDoc.NewElement("ExportDirectoryPath");
    pElement->SetText("C:/Users/Sarah/Desktop/TestRecon");
    pRoot->InsertEndChild(pElement);

    XMLError eResult = xmlDoc.SaveFile(filename.c_str());
    XMLCheckResult(eResult);
    return XML_SUCCESS;
}
int XmlUtils::writeReconstructionTemplate(std::string filename)
{
    XMLDocument xmlDoc;
    XMLNode* pRoot = xmlDoc.NewElement("Root");
    xmlDoc.InsertEndChild(pRoot);

    // First projection 
    XMLElement* pElement = xmlDoc.NewElement("Projection");
    XMLElement* pChild = xmlDoc.NewElement("ImageFilepath");
    pChild->SetText("C:/Users/Sarah/Desktop/TestRecon/image1.ann");
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Second projection ...
/*  pElement = xmlDoc.NewElement("Projection");
    pChild = xmlDoc.NewElement("ImageFilepath");
    pChild->SetText("C:/Users/Sarah/Desktop/TestRecon/image2.ann");
    pElement->InsertEndChild(pChild);

    ... etc.
    */

    // Reconstruction parameters
    pElement = xmlDoc.NewElement("ReconstructionParameters");
    pElement->SetAttribute("type", static_cast<int>(ReconstructionParameters::ReconstructionType::ConstrainedDist));
    pElement->SetAttribute("maxDist", 3.0);
    pElement->SetAttribute("voxelSize", 1.0);
    pChild = xmlDoc.NewElement("VolumeSize");
    pChild->SetAttribute("x", 100);
    pChild->SetAttribute("y", 100);
    pChild->SetAttribute("z", 100);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("VolumeOrigin");
    pChild->SetAttribute("x", 0);
    pChild->SetAttribute("y", 0);
    pChild->SetAttribute("z", 0);
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Reconstruction thresholds (only needed if type is ConstrainedDists)
    pElement = xmlDoc.NewElement("ReconstructionThresholds");
    pElement->SetAttribute("maxDistThreshold", 0.5);
    pElement->SetAttribute("maxRadiusDiff", 0.5);
    pElement->SetAttribute("maxTimeDiff", 10);
    pElement->SetAttribute("maxVelocityComponentDiff", 1.0);
    pRoot->InsertEndChild(pElement);

    // Exporting directory
    pElement = xmlDoc.NewElement("ExportDirectoryPath");
    pElement->SetText("C:/Users/Sarah/Desktop/TestRecon");
    pRoot->InsertEndChild(pElement);

    XMLError eResult = xmlDoc.SaveFile(filename.c_str());
    XMLCheckResult(eResult);
    return XML_SUCCESS;
}
int XmlUtils::writeValidationTemplate(std::string filename)
{
    XMLDocument xmlDoc;
    XMLNode* pRoot = xmlDoc.NewElement("Root");
    xmlDoc.InsertEndChild(pRoot);

    XMLElement* pElement = xmlDoc.NewElement("BravaFilepath");
    pElement->SetText("C:/Users/Sarah/Desktop/TestRecon/BG01_ColorCoded.swc");
    pRoot->InsertEndChild(pElement);

    pElement = xmlDoc.NewElement("MaxImageDimension");
    pElement->SetAttribute("maxDim", 256);
    pRoot->InsertEndChild(pElement);

    pElement = xmlDoc.NewElement("MaxDistanceInImage");
    pElement->SetAttribute("maxDist", 3.0);
    pRoot->InsertEndChild(pElement);

    pElement = xmlDoc.NewElement("Circulation");
    pElement->SetAttribute("id", 2);
    pRoot->InsertEndChild(pElement);
    pElement = xmlDoc.NewElement("Circulation");
    pElement->SetAttribute("id", 3);
    pRoot->InsertEndChild(pElement);
    pElement = xmlDoc.NewElement("Circulation");
    pElement->SetAttribute("id", 6);
    pRoot->InsertEndChild(pElement);

    // First projection (AP)
    pElement = xmlDoc.NewElement("CameraData");
    XMLElement* pChild = xmlDoc.NewElement("CameraFocalLength");
    pChild->SetAttribute("f", 1000);    // In mm
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("DirectionToCamera");
    pChild->SetAttribute("x", 0.0);
    pChild->SetAttribute("y", -1.0);    // AP view
    pChild->SetAttribute("z", 0.0);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("CameraUpVector");
    pChild->SetAttribute("x", 0.0);
    pChild->SetAttribute("y", 0.0);
    pChild->SetAttribute("z", 1.0);     // IS direction
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Second projection (RL)
    pElement = xmlDoc.NewElement("CameraData");
    pChild = xmlDoc.NewElement("CameraFocalLength");
    pChild->SetAttribute("f", 1000);    // In mm
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("DirectionToCamera");
    pChild->SetAttribute("x", 1.0);
    pChild->SetAttribute("y", 0.0);    // RL view
    pChild->SetAttribute("z", 0.0);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("CameraUpVector");
    pChild->SetAttribute("x", 0.0);
    pChild->SetAttribute("y", 0.0);
    pChild->SetAttribute("z", 1.0);     // IS direction
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Reconstruction threshold ranges
    pElement = xmlDoc.NewElement("ReconstructionThresholdRanges");
    pChild = xmlDoc.NewElement("DistRange");        // mm from centerline
    pChild->SetAttribute("min", 1);
    pChild->SetAttribute("max", 2);
    pChild->SetAttribute("numSamples", 1);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("RadiusRange");      // in mm
    pChild->SetAttribute("min", 5);
    pChild->SetAttribute("max", 10);
    pChild->SetAttribute("numSamples", 1);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("TimeRange");        // [0, maxTime]
    pChild->SetAttribute("min", 5);
    pChild->SetAttribute("max", 10);
    pChild->SetAttribute("numSamples", 1);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("SpeedRange");       // [0, 1] for Brava data   
    pChild->SetAttribute("min", 0.5);
    pChild->SetAttribute("max", 1);
    pChild->SetAttribute("numSamples", 1);
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Second projection ...
    //pElement = xmlDoc.NewElement("ReconstructionThresholdRanges");
    //pChild = xmlDoc.NewElement("DistRange");        // mm from centerline
    //pChild->SetAttribute("min", 1);
    //pChild->SetAttribute("max", 2);
    //pChild->SetAttribute("numSamples", 1);
    //pElement->InsertEndChild(pChild);
    //... etc.
    
    // Exporting directory
    pElement = xmlDoc.NewElement("ExportFilePath");
    pElement->SetText("C:/Users/Sarah/Desktop/TestRecon/metrics.txt");
    pRoot->InsertEndChild(pElement);

    XMLError eResult = xmlDoc.SaveFile(filename.c_str());
    XMLCheckResult(eResult);
    return XML_SUCCESS;
}

int XmlUtils::writeReconstruction(std::string filename, std::vector<std::string> imageFilenames,
    ReconstructionParameters reconParams, ReconstructionThresholds reconThresholds)
{
    XMLDocument xmlDoc;
    XMLNode* pRoot = xmlDoc.NewElement("Root");
    xmlDoc.InsertEndChild(pRoot);

    for (int i = 0; i < imageFilenames.size(); i++) {
        XMLElement* pElement = xmlDoc.NewElement("Projection");
        XMLElement* pChild = xmlDoc.NewElement("ImageFilepath");
        pChild->SetText(imageFilenames[i].c_str());
        pElement->InsertEndChild(pChild);
        pRoot->InsertEndChild(pElement);
    }

    // Reconstruction parameters
    XMLElement* pElement = xmlDoc.NewElement("ReconstructionParameters");
    pElement->SetAttribute("type", static_cast<int>(reconParams.type()));
    pElement->SetAttribute("maxDist", reconParams.maxDist());
    pElement->SetAttribute("voxelSize", reconParams.voxelSize());
    XMLElement* pChild = xmlDoc.NewElement("VolumeSize");
    pChild->SetAttribute("x", reconParams.volumeSize()[0]);
    pChild->SetAttribute("y", reconParams.volumeSize()[1]);
    pChild->SetAttribute("z", reconParams.volumeSize()[2]);
    pElement->InsertEndChild(pChild);
    pChild = xmlDoc.NewElement("VolumeOrigin");
    pChild->SetAttribute("x", reconParams.volumeOrigin().x());
    pChild->SetAttribute("y", reconParams.volumeOrigin().y());
    pChild->SetAttribute("z", reconParams.volumeOrigin().z());
    pElement->InsertEndChild(pChild);
    pRoot->InsertEndChild(pElement);

    // Reconstruction thresholds
    pElement = xmlDoc.NewElement("ReconstructionThresholds");
    pElement->SetAttribute("maxDistThreshold", reconThresholds.maxDist());
    pElement->SetAttribute("maxRadiusDiff", reconThresholds.maxRadiusDiff());
    pElement->SetAttribute("maxTimeDiff", reconThresholds.maxTimeDiff());
    pElement->SetAttribute("maxVelocityComponentDiff", reconThresholds.maxSharedVelCompDiff());
    pRoot->InsertEndChild(pElement);

    // Exporting directory
    pElement = xmlDoc.NewElement("ExportDirectoryPath");
    pElement->SetText("C:/Users/Sarah/Desktop/TestRecon");
    pRoot->InsertEndChild(pElement);

    XMLError eResult = xmlDoc.SaveFile(filename.c_str());
    XMLCheckResult(eResult);
    return XML_SUCCESS;
}

int XmlUtils::readSyntheticGenerationXml(std::string filename)
{
    m_syntheticGenerationRead = false;

    XMLDocument xmlDoc;
    int eResult = xmlDoc.LoadFile(filename.c_str());
    XMLCheckResult(eResult);

    // Parse xml file
    XMLNode* pRoot = xmlDoc.FirstChild();
    if (pRoot == nullptr) return XML_ERROR_FILE_READ_ERROR;

    XMLElement* pElement = pRoot->FirstChildElement("BravaFilepath");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    m_modelFilepath = pElement->GetText();

    pElement = pRoot->FirstChildElement("MaxImageDimension");
    pElement->QueryIntAttribute("maxDim", &m_maxDim);
    pRoot->InsertEndChild(pElement);

    pElement = pRoot->FirstChildElement("MaxDistanceInImage");
    pElement->QueryFloatAttribute("maxDist", &m_maxDistInImage);
    pRoot->InsertEndChild(pElement);

    pElement = pRoot->FirstChildElement("Circulation");
    m_circulations.clear();
    while (pElement) {
        int id;
        pElement->QueryIntAttribute("id", &id);
        m_circulations.push_back(id);
        pElement = pElement->NextSiblingElement("Circulation");
    }

    // Parse the elements containing projection data
    pElement = pRoot->FirstChildElement("CameraData");
    m_cameraData.clear();
    while (pElement) {
        CameraData cd;
        XMLElement* pChild = pElement->FirstChildElement("CameraFocalLength");
        pChild->QueryFloatAttribute("f", &(cd.focalLength));
        pChild = pElement->FirstChildElement("DirectionToCamera");
        pChild->QueryFloatAttribute("x", &(cd.cameraDir[0]));
        pChild->QueryFloatAttribute("y", &(cd.cameraDir[1]));
        pChild->QueryFloatAttribute("z", &(cd.cameraDir[2]));
        pChild = pElement->FirstChildElement("CameraUpVector");
        pChild->QueryFloatAttribute("x", &(cd.cameraUp[0]));
        pChild->QueryFloatAttribute("y", &(cd.cameraUp[1]));
        pChild->QueryFloatAttribute("z", &(cd.cameraUp[2]));
        m_cameraData.push_back(cd);
        pElement = pElement->NextSiblingElement("CameraData");
    }

    // Export options
    int exportVolBool(0);
    int exportImagesBool(0);
    pElement = pRoot->FirstChildElement("ExportOptions");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    pElement->QueryIntAttribute("doExportVolumes", &exportVolBool);
    pElement->QueryIntAttribute("doExportAllImages", &exportImagesBool);
    m_doExportVols = (exportVolBool == 0) ? false : true;
    m_doExportAllImages = (exportImagesBool == 0) ? false : true;

    // Export directory
    pElement = pRoot->FirstChildElement("ExportDirectoryPath");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    m_genExportPath = pElement->GetText();

    m_syntheticGenerationRead = true;
    return XML_SUCCESS;
}
int XmlUtils::readReconstructionXml(std::string filename)
{
    m_reconstructionRead = true;

    XMLDocument xmlDoc;
    int eResult = xmlDoc.LoadFile(filename.c_str());
    XMLCheckResult(eResult);

    // Parse xml file
    XMLNode* pRoot = xmlDoc.FirstChild();
    if (pRoot == nullptr) return XML_ERROR_FILE_READ_ERROR;

    // Parse the elements containing projection data
    XMLElement* pElement = pRoot->FirstChildElement("Projection");
    m_reconImageFilenames.clear();
    while (pElement) {
        XMLElement* pChild = pElement->FirstChildElement("ImageFilepath");
        m_reconImageFilenames.push_back(pChild->GetText());
        pElement = pElement->NextSiblingElement("Projection");
    }

    int type = 0;
    float maxDist = 10;
    float voxelSize = 1;
    std::array<int, 3> volumeSize = { 0, 0, 0 };
    float volumeOrigin[3] = { 0, 0, 0 };
    pElement = pRoot->FirstChildElement("ReconstructionParameters");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    pElement->QueryIntAttribute("type", &type);
    pElement->QueryFloatAttribute("maxDist", &maxDist);
    pElement->QueryFloatAttribute("voxelSize", &voxelSize);
    XMLElement* pChild = pElement->FirstChildElement("VolumeSize");
    if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
    pChild->QueryIntAttribute("x", &(volumeSize[0]));
    pChild->QueryIntAttribute("y", &(volumeSize[1]));
    pChild->QueryIntAttribute("z", &(volumeSize[2]));
    pChild = pElement->FirstChildElement("VolumeOrigin");
    if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
    pChild->QueryFloatAttribute("x", &(volumeOrigin[0]));
    pChild->QueryFloatAttribute("y", &(volumeOrigin[1]));
    pChild->QueryFloatAttribute("z", &(volumeOrigin[2]));
    if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
    ReconstructionParameters::ReconstructionType reconType = static_cast<ReconstructionParameters::ReconstructionType>(type);
    QVector3D origin(volumeOrigin[0], volumeOrigin[1], volumeOrigin[2]);
    ReconstructionParameters reconParams(reconType, maxDist, voxelSize, volumeSize, origin);
    m_reconParams = reconParams;

    // Reconstruction thresholds
    float maxDistThreshold = 0;
    float maxRadDiff = 0;
    float maxTimeDiff = 0;
    float maxSpeedDiff = 0;
    pElement = pRoot->FirstChildElement("ReconstructionThresholds");
    pElement->QueryFloatAttribute("maxDistThreshold", &maxDistThreshold);
    pElement->QueryFloatAttribute("maxRadiusDiff", &maxRadDiff);
    pElement->QueryFloatAttribute("maxTimeDiff", &maxTimeDiff);
    pElement->QueryFloatAttribute("maxVelocityComponentDiff", &maxSpeedDiff);
    ReconstructionThresholds rt(maxDistThreshold, maxRadDiff, maxTimeDiff, maxSpeedDiff);
    m_reconThresholds = rt;

    pElement = pRoot->FirstChildElement("ExportDirectoryPath");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    m_reconExportPath = pElement->GetText();

    m_reconstructionRead = true;
    return XML_SUCCESS;
}

int XmlUtils::readValidationXml(std::string filename)
{
    m_validationRead = false;

    XMLDocument xmlDoc;
    int eResult = xmlDoc.LoadFile(filename.c_str());
    XMLCheckResult(eResult);

    // Parse xml file
    XMLNode* pRoot = xmlDoc.FirstChild();
    if (pRoot == nullptr) return XML_ERROR_FILE_READ_ERROR;

    XMLElement* pElement = pRoot->FirstChildElement("BravaFilepath");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    m_modelFilepath = pElement->GetText();

    pElement = pRoot->FirstChildElement("MaxImageDimension");
    pElement->QueryIntAttribute("maxDim", &m_maxDim);
    pRoot->InsertEndChild(pElement);

    pElement = pRoot->FirstChildElement("MaxDistanceInImage");
    pElement->QueryFloatAttribute("maxDist", &m_maxDistInImage);
    pRoot->InsertEndChild(pElement);

    pElement = pRoot->FirstChildElement("Circulation");
    m_circulations.clear();
    while (pElement) {
        int id;
        pElement->QueryIntAttribute("id", &id);
        m_circulations.push_back(id);
        pElement = pElement->NextSiblingElement("Circulation");
    }

    // Parse the elements containing projection data
    pElement = pRoot->FirstChildElement("CameraData");
    m_cameraData.clear();
    while (pElement) {
        CameraData cd;
        XMLElement* pChild = pElement->FirstChildElement("CameraFocalLength");
        pChild->QueryFloatAttribute("f", &(cd.focalLength));
        pChild = pElement->FirstChildElement("DirectionToCamera");
        pChild->QueryFloatAttribute("x", &(cd.cameraDir[0]));
        pChild->QueryFloatAttribute("y", &(cd.cameraDir[1]));
        pChild->QueryFloatAttribute("z", &(cd.cameraDir[2]));
        pChild = pElement->FirstChildElement("CameraUpVector");
        pChild->QueryFloatAttribute("x", &(cd.cameraUp[0]));
        pChild->QueryFloatAttribute("y", &(cd.cameraUp[1]));
        pChild->QueryFloatAttribute("z", &(cd.cameraUp[2]));
        m_cameraData.push_back(cd);
        pElement = pElement->NextSiblingElement("CameraData");
    }

    // Parse the elements containing threshold ranges
    // Reconstruction threshold ranges
    pElement = pRoot->FirstChildElement("ReconstructionThresholdRanges");
    m_validationThresholds.clear();
    while (pElement) {
        ThresholdRange distRange;
        ThresholdRange radiusRange;
        ThresholdRange timeRange;
        ThresholdRange speedRange;
        if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
        XMLElement* pChild = pElement->FirstChildElement("DistRange");
        if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
        pChild->QueryFloatAttribute("min", &(distRange.min));
        pChild->QueryFloatAttribute("max", &distRange.max);
        pChild->QueryIntAttribute("numSamples", &distRange.numSamples);
        pChild = pElement->FirstChildElement("RadiusRange");
        if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
        pChild->QueryFloatAttribute("min", &radiusRange.min);
        pChild->QueryFloatAttribute("max", &radiusRange.max);
        pChild->QueryIntAttribute("numSamples", &radiusRange.numSamples);
        pChild = pElement->FirstChildElement("TimeRange");
        if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
        pChild->QueryFloatAttribute("min", &timeRange.min);
        pChild->QueryFloatAttribute("max", &timeRange.max);
        pChild->QueryIntAttribute("numSamples", &timeRange.numSamples);
        pChild = pElement->FirstChildElement("SpeedRange");
        if (pChild == nullptr) return XML_ERROR_PARSING_ELEMENT;
        pChild->QueryFloatAttribute("min", &speedRange.min);
        pChild->QueryFloatAttribute("max", &speedRange.max);
        pChild->QueryIntAttribute("numSamples", &speedRange.numSamples);
        std::vector<ReconstructionThresholds> vt = thresholdsFromRanges(distRange, radiusRange, timeRange, speedRange);
        m_validationThresholds.insert(m_validationThresholds.end(), vt.begin(), vt.end());
        pElement = pElement->NextSiblingElement("ReconstructionThresholdRanges");
    }

    pElement = pRoot->FirstChildElement("ExportFilePath");
    if (pElement == nullptr) return XML_ERROR_PARSING_ELEMENT;
    m_valExportPath = pElement->GetText();

    // Set booleans so data can be accessed
    m_syntheticGenerationRead = true;
    m_reconstructionRead = true;
    m_validationRead = true;
    return XML_SUCCESS;
}

// Access state read from a synthetic generation xml file
std::string XmlUtils::modelFilename()
{
    std::string filename = "";
    if (m_syntheticGenerationRead) filename = m_modelFilepath;
    return filename;
}
int XmlUtils::maxDim()
{
    int maxDim = 0;
    if (m_syntheticGenerationRead) maxDim = m_maxDim;
    return maxDim;
}
float XmlUtils::maxDistInImage()
{
    float maxDist = 0;
    if (m_syntheticGenerationRead) maxDist = m_maxDistInImage;
    return maxDist;
}
bool XmlUtils::doExportVolumes()
{
    bool doExport(false);
    if (m_syntheticGenerationRead) doExport = m_doExportVols;
    return doExport;
}
bool XmlUtils::doExportAllImages()
{
    bool doExport(false);
    if (m_syntheticGenerationRead) doExport = m_doExportAllImages;
    return doExport;
}
std::vector<int> XmlUtils::circulations()
{
    std::vector<int> circulations;
    if (m_syntheticGenerationRead) circulations = m_circulations;
    return circulations;
}
std::vector<Projector::CameraParameters> XmlUtils::genCameraParams()
{
    std::vector<Projector::CameraParameters> cameraParams;
    if (m_syntheticGenerationRead) {
        for (std::vector<CameraData>::iterator it = m_cameraData.begin(); it < m_cameraData.end(); it++) {
            // Compute camera parameters assuming the camera looks through the world
            // origin and focuses at a point equidistant from the origin.
            Projector::CameraParameters cp;
            cp.focalLength = it->focalLength;
            QVector3D origin(0, 0, 0);
            QVector3D dirToCamera = (it->cameraDir).normalized();
            QVector3D cameraPos = origin + dirToCamera * 0.5 * it->focalLength;
            cp.worldToCameraMatrix = worldToCameraMatrix(cameraPos, origin, it->cameraUp);
            cameraParams.push_back(cp);
        }
    }
    return cameraParams;
}
std::string XmlUtils::genExportDir()
{
    return m_genExportPath;
}

// Access state read from a reconstruction xml file
std::vector<std::string> XmlUtils::imageFilenames()
{
    std::vector<std::string> imageFilenames;
    if (m_reconstructionRead) imageFilenames = m_reconImageFilenames;
    return imageFilenames;
}
ReconstructionParameters XmlUtils::reconParams()
{
    ReconstructionParameters params;
    if (m_reconstructionRead) params = m_reconParams;
    return params;
}
ReconstructionThresholds XmlUtils::reconThresholds()
{
    ReconstructionThresholds thresholds;
    if (m_reconstructionRead) thresholds = m_reconThresholds;
    return thresholds;
}
std::string XmlUtils::reconExportDir()
{
    return m_reconExportPath;
}
std::vector<ReconstructionThresholds> XmlUtils::validationThresholds()
{
    return m_validationThresholds;
}
std::string XmlUtils::validateExportFilename()
{
    return m_valExportPath;
}


//
// Private
//
std::vector<ReconstructionThresholds> XmlUtils::thresholdsFromRanges(ThresholdRange dist,
    ThresholdRange radius, ThresholdRange time, ThresholdRange speed)
{
    if (dist.max <= dist.min || dist.numSamples < 1) dist.numSamples = 1;
    float distIncr = (dist.numSamples < 2) ? 0 : (dist.max - dist.min) / (dist.numSamples - 1);
    if (radius.max <= radius.min || radius.numSamples < 1) radius.numSamples = 1;
    float radiusIncr = (radius.numSamples < 2) ? 0 : (radius.max - radius.min) / (radius.numSamples - 1);
    if (time.max <= time.min || time.numSamples < 1) time.numSamples = 1;
    float timeIncr = (time.numSamples < 2) ? 0 : (time.max - time.min) / (time.numSamples - 1);
    if (speed.max <= speed.min || speed.numSamples < 1) speed.numSamples = 1;
    float speedIncr = (speed.numSamples < 2) ? 0 : (speed.max - speed.min) / (speed.numSamples - 1);

    std::vector<ReconstructionThresholds> t;
    for (int distIdx = 0; distIdx < dist.numSamples; distIdx++) {
        float distT = dist.min + distIdx * distIncr;
        for (int radiusIdx = 0; radiusIdx < radius.numSamples; radiusIdx++) {
            float radiusT = radius.min + radiusIdx * radiusIncr;
            for (int timeIdx = 0; timeIdx < time.numSamples; timeIdx++) {
                float timeT = time.min + timeIdx * timeIncr;
                for (int speedIdx = 0; speedIdx < speed.numSamples; speedIdx++) {
                    float speedT = speed.min + speedIdx * speedIncr;
                    t.push_back(ReconstructionThresholds(distT, radiusT, timeT, speedT));
                }
            }
        }
    }
    return t;
}

QMatrix4x4 XmlUtils::worldToCameraMatrix(const QVector3D& cameraPos, const QVector3D& targetPos,
    const QVector3D& desiredUp)
{
    QVector3D targetToCamera = cameraPos - targetPos;
    targetToCamera.normalize();
    QVector3D left = QVector3D::crossProduct(desiredUp, targetToCamera);
    left.normalize();
    QVector3D up = QVector3D::crossProduct(targetToCamera, left);

    QMatrix4x4 R; // Rotate target to camera view, initialized to identity
    R(0, 0) = left.x();
    R(1, 0) = left.y();
    R(2, 0) = left.z();
    R(0, 1) = up.x();
    R(1, 1) = up.y();
    R(2, 1) = up.z();
    R(0, 2) = targetToCamera.x();
    R(1, 2) = targetToCamera.y();
    R(2, 2) = targetToCamera.z();

    QMatrix4x4 lookat;
    lookat = R.transposed();
    lookat.translate(-cameraPos);
    return lookat; // World to camera transform
}

void XmlUtils::removeDoubleBackslash(std::string& string)
{
    std::string doubleBackslash("\\\\");
    std::string backslash("\\");
    size_t pos = 0;
    while ((pos = string.find(doubleBackslash, pos)) != std::string::npos) {
        string.replace(pos, doubleBackslash.length(), backslash);
        pos += backslash.length();
    }
}
