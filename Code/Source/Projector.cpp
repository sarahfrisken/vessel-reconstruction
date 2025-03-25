//
// Projector.cpp
// Implementation of Projector.
//

#include "Projector.h"

// 
// Public
//
Projector::Projector(Projector::CameraParameters cameraParams,
	Projector::DetectorParameters detectorParams) :
	m_cameraParams(cameraParams),
	m_detectorParams(detectorParams)
{
	m_worldToCamera = cameraParams.worldToCameraMatrix;

	computeProjectionTransform();
	computePixelTransform();
}
Projector::~Projector()
{
}

QVector4D Projector::worldToImage(QVector3D worldPoint)
{
	// Transform the point from world to camera coordinates and project it to the image 
	QVector4D hWorldPoint = { worldPoint, 1 };
	QVector4D imagePoint = m_cameraToImage * m_worldToCamera * hWorldPoint;

	// After the perspective divide, (x, y) is the point in image coordinates, z = f, 
	// the focal distance (i.e., the distance from the camera to the image), and w = 1.
	imagePoint /= imagePoint[3];
	return imagePoint;
}
QVector2D Projector::worldToPixel(QVector3D worldPoint)
{
	QVector4D imagePoint = worldToImage(worldPoint);
	return QVector2D(m_imageToPixel * imagePoint);
}
float Projector::worldLengthToPixelScale(QVector3D worldPoint)
{
	// Project a small vector at the world point that is perpendicular to the camera 
	// direction onto the image and return the scale of that vector's length
	float smallOffset = 0.1;
	QVector4D cameraPoint = worldToCamera(worldPoint);
	QVector4D offsetCameraPoint = cameraPoint + QVector4D(smallOffset, smallOffset, 0, 0);

	QVector2D pixelPoint = cameraToPixel(cameraPoint);
	QVector2D offsetPixelPoint = cameraToPixel(offsetCameraPoint);
	float lengthScale = (offsetPixelPoint - pixelPoint).length() / smallOffset;
	return lengthScale;
}
float Projector::pixelLengthToWorldScale(QVector3D worldPoint)
{
	float scale = 1.0f / worldLengthToPixelScale(worldPoint);
	return scale;
}
QVector2D Projector::worldVectorToPixel(QVector3D worldPoint, QVector3D vector)
{
	// Project the vector endpoints onto the image to determine the projected vector
	QVector2D pixelStart = worldToPixel(worldPoint);
	QVector2D pixelEnd = worldToPixel(worldPoint + vector);
	return (pixelEnd - pixelStart);
}

// For testing
QVector4D Projector::cameraToWorld(QVector3D worldPoint)
{
	QMatrix4x4 cameraToWorld = m_worldToCamera.inverted();
	QVector4D hWorldPoint = { worldPoint, 1 };
	return cameraToWorld * hWorldPoint;
}

// 
// Private
//
void Projector::computeProjectionTransform()
{
	m_cameraToImage.setToIdentity();
	m_cameraToImage.scale(m_cameraParams.focalLength);

	// Perspective projection
	m_cameraToImage(3,2) = 1;
	m_cameraToImage(3,3) = 0;
}
void Projector::computePixelTransform()
{
	assert(m_detectorParams.pixelSize.x() != 0);
	assert(m_detectorParams.pixelSize.y() != 0);

	m_imageToPixel.setToIdentity();
	QVector3D originOffset = { (float)m_detectorParams.size.width() * 0.5f,
		 (float)m_detectorParams.size.height() * 0.5f, 0 };
	m_imageToPixel.translate(originOffset);
	m_imageToPixel.scale(1.0f / m_detectorParams.pixelSize.x(),
		1.0f / m_detectorParams.pixelSize.y(), 1);
}

QVector4D Projector::worldToCamera(QVector3D worldPoint)
{
	QVector4D hWorldPoint = { worldPoint, 1 };
	return m_worldToCamera * hWorldPoint;
}
QVector4D Projector::cameraToImage(QVector4D cameraPoint)
{
	return m_cameraToImage * cameraPoint;
}
QVector2D Projector::cameraToPixel(QVector4D cameraPoint)
{
	QVector4D imagePoint = cameraToImage(cameraPoint);

	// After the perspective divide, (x, y) is the point in image coordinates, z = f, 
	// the focal distance (i.e., the distance from the camera to the image), and w = 1.
	imagePoint /= imagePoint[3];
	return QVector2D(m_imageToPixel * imagePoint);
}