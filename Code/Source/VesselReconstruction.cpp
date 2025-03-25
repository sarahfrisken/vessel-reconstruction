//
// VesselReconstruction.cpp
// Implementation of VesselReconstruction.
//

#include "VesselReconstruction.h"

#include "XmlUtils.h"
#include "BravaUtils.h"
#include "Projector.h"
#include "Reconstructor.h"
#include "Validator.h"

#include <QtWidgets>

#include <assert.h>
#include <iostream>
#include <fstream>

VesselReconstruction::VesselReconstruction(QWidget *parent) :
    QMainWindow(parent),
    m_voxelSize(0.2),
    m_maxDist(3.0)
{
    ui.setupUi(this);

    connect(ui.generateButton, &QPushButton::released, this, &VesselReconstruction::onGenerateXmlTemplates);
    connect(ui.preprocessButton, &QPushButton::released, this, &VesselReconstruction::onPreprocessBravaModels);
    connect(ui.reconstructButton, &QPushButton::released, this, &VesselReconstruction::onReconstruct);
    connect(ui.validateButton, &QPushButton::released, this, &VesselReconstruction::onValidate);
}

VesselReconstruction::~VesselReconstruction()
{
}

void VesselReconstruction::onGenerateXmlTemplates()
{
    QString filedir = QFileDialog::getExistingDirectory(this, tr("Directory for storing templates"), "/home", 
        QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    if (!filedir.isEmpty() && !filedir.isNull()) {
        XmlUtils::writeGenerateFromBravaTemplate((filedir + "/preprocess.xml").toStdString());
        XmlUtils::writeReconstructionTemplate((filedir + "/reconstruct.xml").toStdString());
        XmlUtils::writeValidationTemplate((filedir + "/validate.xml").toStdString());
    }
}

void VesselReconstruction::onPreprocessBravaModels()
{
    QString filename = QFileDialog::getOpenFileName(0, tr("XML file controlling generation"), "/home", tr("*.xml"));
    XmlUtils xmlReader;
    int error = xmlReader.readSyntheticGenerationXml(filename.toStdString().c_str());
    if (error) {
        std::cout << "Bad generation file" << "\n";
        return;
    }

    // Get Brava file and data
    BravaUtils brava;
    brava.importBravaData(xmlReader.modelFilename().c_str());

    // Get generation data. Brava circulations: (2, 3, 6) for right ICA, (4, 5, 7) for left ICA.
    std::vector<int> circulations = xmlReader.circulations();
    std::vector<VesselModel::Line> lines = brava.vesselLines(circulations);
    VesselModel model(lines);

    // Export path
    std::string dir = xmlReader.genExportDir();

    // Export the model
    if (xmlReader.doExportVolumes()) {
        model.exportSignedDistToVesselVolume((dir + "/distToVessel").c_str(), m_voxelSize, m_maxDist);
        model.exportDistToCenterlineVolume((dir + "/distToCenterline").c_str(), m_voxelSize, m_maxDist);
    }

    // Create and export the specified annotated projection images for this model
    int maxDim = xmlReader.maxDim();
    float maxDist = xmlReader.maxDistInImage();
    std::vector<Projector::CameraParameters> cameraParams = xmlReader.genCameraParams();
    std::vector<std::string> filenames;
    for (int i = 0; i < cameraParams.size(); i++) {
        // Determine the optimal image parameters and perform the projection
        Projector::DetectorParameters optImageParms = model.getOptimalImageParams(cameraParams[i], maxDim, maxDist);
        std::shared_ptr<ProjectionImage> projImage = std::make_shared<ProjectionImage>(cameraParams[i], optImageParms);
        model.renderProjectionImage(projImage, maxDist);

        // Export projection images. Uses JPEG images because QImage save for png and bmp don't
        // doesn't seem to work with alpha, which is used to encode distances.
        std::string prefix = ("/view" + std::to_string(i) + "_");
        if (xmlReader.doExportAllImages()) {
            projImage->exportProjImages((dir + prefix).c_str());
        }
        std::string annImageFilename = dir + prefix + "annotation.ann";
        projImage->exportAnnotation(annImageFilename.c_str());

        // Prepare to export reconstruction xml file
        filenames.push_back(annImageFilename);
    }

     // Create the reconstruction XML file for this set of projections
    float voxelSize = m_voxelSize;
    VesselModel::BBox bbox = model.boundingBox();
    QVector3D origin = bbox.min;
    QVector3D fSize = (bbox.max - bbox.min) / voxelSize;
    std::array<int, 3> size = { int(fSize.x() + 0.5), int(fSize.y() + 0.5), int(fSize.z() + 0.5) };
    ReconstructionParameters::ReconstructionType reconType = ReconstructionParameters::ReconstructionType::ConstrainedDist;
    ReconstructionParameters reconParams(reconType, m_maxDist, voxelSize, size, origin);

   // Default type and thresholds
    float maxDistThreshold = m_voxelSize;
    float maxRadiusDiff = 0.5;
    float maxTimeDiff = 10;
    float maxSpeedDiff = 1;
    ReconstructionThresholds reconThresholds(maxDistThreshold, maxRadiusDiff, maxTimeDiff, maxSpeedDiff);

    XmlUtils::writeReconstruction(dir + "/reconstruct.xml", filenames, reconParams, reconThresholds);
}
void VesselReconstruction::onReconstruct()
{
    QString filename = QFileDialog::getOpenFileName(0, tr("XML file controlling reconstruction"), "reconstruct.xml", tr("*.xml"));
    XmlUtils xmlReader;
    int error = xmlReader.readReconstructionXml(filename.toStdString().c_str());
    if (error) {
        std::cout << "Bad reconstruction file" << "\n";
        return;
    }

    // Get projection images
    std::vector<std::string> filenames = xmlReader.imageFilenames();
    std::vector<std::shared_ptr<ProjectionImage>> projImages;
    for (std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++) {
        projImages.emplace_back(std::move(std::make_shared<ProjectionImage>(it->c_str())));
    }
    
    // Get reconstruction parameters
    ReconstructionParameters reconParams = xmlReader.reconParams();
    ReconstructionThresholds reconThresholds = xmlReader.reconThresholds();

    // Use test thresolds and perform the reconstruction
    Reconstructor reconstructor(projImages, reconParams);
    reconstructor.reconstructAndExport((xmlReader.reconExportDir() + "/vesselnessVol").c_str(), 
        reconThresholds);
}
void VesselReconstruction::onValidate()
{
    QString filename = QFileDialog::getOpenFileName(0, tr("XML file controlling validation"), "validate.xml", tr("*.xml"));
    XmlUtils xmlReader;
    int error = xmlReader.readValidationXml(filename.toStdString().c_str());

    // Generate model distance field and projection images
    BravaUtils brava;
    brava.importBravaData(xmlReader.modelFilename().c_str());
    std::vector<int> circulations = xmlReader.circulations();
    std::vector<VesselModel::Line> lines;
    lines = brava.vesselLines(circulations);
    VesselModel model(lines);

    // Generate projection images
    int maxDim = xmlReader.maxDim();
    float maxDist = xmlReader.maxDistInImage();
    std::vector<Projector::CameraParameters> cameraParams = xmlReader.genCameraParams();
    if (cameraParams.size() < 2) return;
    std::vector<std::shared_ptr<ProjectionImage>> projImages;
    for (int i = 0; i < cameraParams.size(); i++) {
        // Determine the optimal image parameters and perform the projection
        Projector::DetectorParameters optImageParms = model.getOptimalImageParams(cameraParams[i], maxDim, maxDist);
        projImages.emplace_back(std::move(std::make_shared<ProjectionImage>(cameraParams[i], optImageParms)));
        model.renderProjectionImage(projImages.back(), maxDist);
    }

    // Get validation data
    ReconstructionParameters::ReconstructionType type;
    std::vector<ReconstructionThresholds> thresholds;
    if (cameraParams.size() > 2) {
        type = ReconstructionParameters::ReconstructionType::MaxDist;
    }
    else {
        type = ReconstructionParameters::ReconstructionType::ConstrainedDist;
    }
    thresholds = xmlReader.validationThresholds();

    // Set reconstruction parameters
    float voxelSize = m_voxelSize;
    VesselModel::BBox bbox = model.boundingBox();
    QVector3D origin = bbox.min;
    QVector3D fSize = (bbox.max - bbox.min) / voxelSize;
    std::array<int, 3> size = { int(fSize.x() + 0.5), int(fSize.y() + 0.5), int(fSize.z() + 0.5) };
    ReconstructionParameters reconParams(type, m_maxDist, voxelSize, size, origin);

    // Perform the validation and export the results
        // Create the output stream 
    std::ofstream fstream(xmlReader.validateExportFilename().c_str());
    try {
        if (!fstream) {
            throw std::runtime_error("Cannot create the metrics export.");
        }
    }
    catch (const std::exception& e) {
        std::cout << "Exception " << e.what() << std::endl;
        return;
    }

    Reconstructor reconstructor(projImages, reconParams);
    std::shared_ptr<Volume> distVol = model.distToCenterlineVolume();
    std::vector<ValidationMetrics> metrics = reconstructor.reconstructAndValidate(distVol, thresholds);
    for (int i = 0; i < thresholds.size(); i++) {
        fstream << thresholds[i].maxDist() << ", "
            << thresholds[i].maxRadiusDiff() << ", "
            << thresholds[i].maxTimeDiff() << ", "
            << thresholds[i].maxSharedVelCompDiff() << ", "
            << metrics[i].truePos() << ", "
            << metrics[i].trueNeg() << ", "
            << metrics[i].falsePos() << ", "
            << metrics[i].falseNeg() << "\n";
    }
}


