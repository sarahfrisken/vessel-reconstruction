//
// Projector.h
// Utility for simulating the projection of a 3D vessel model onto a 2D image.
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

#include <QSize>
#include <QVector2D>
#include <QVector3D>
#include <QMatrix4x4>

class Projector
{
public:
	struct CameraParameters {
		float focalLength;
		QMatrix4x4 worldToCameraMatrix;
	};
	struct DetectorParameters {
		QSize size;
		QVector2D pixelSize;
	};

	Projector(CameraParameters cameraParams, DetectorParameters detectorParams);
	~Projector();

	QVector4D worldToImage(QVector3D worldPoint);
	QVector2D worldToPixel(QVector3D worldPoint);
	float worldLengthToPixelScale(QVector3D worldPoint);					// Lengths projected onto the image are position dependent
	float pixelLengthToWorldScale(QVector3D worldPoint);
	QVector2D worldVectorToPixel(QVector3D worldPoint, QVector3D vector);   // Projects a 3D vector onto the image

	// For testing
	QVector4D cameraToWorld(QVector3D worldPoint);

private:
	CameraParameters m_cameraParams;
	DetectorParameters m_detectorParams;

	// 4x4 matrices are not typically used for these transforms. They are used here
	// because Qt supports matrix-vector multiplication for 4x4 matricies. 
	QMatrix4x4 m_worldToCamera;
	QMatrix4x4 m_cameraToImage;
	QMatrix4x4 m_imageToPixel;
	void computeProjectionTransform();
	void computePixelTransform();

	QVector4D worldToCamera(QVector3D worldPoint);
	QVector4D cameraToImage(QVector4D cameraPoint);
	QVector2D cameraToPixel(QVector4D cameraPoint);

	//    The transform from a voxel (i, j, k) to world coordinates is given by
	//
	//							|vx   0    0    ox|
	//    Mat3x3 voxelToWorld =	|0    vy   0    oy|,
	//							|0    0    vz   oz|
	//							|0    0    0     1|
	//    where (vx, vy, vz) = voxelSize, (ox, oy, oz) = volume origin in world 
	//    coordinates, volume is assumed aligned with world coordinates. 
	//    If volume size is (Nx, Ny, Nz) and volume is centered at the world origin 
	//    then (ox, oy, oz) = (-vx * Nx/2, -vy * Ny/2, -vz * Nz/2).
	//
	//    The canonical camera looks at the origin facing in the positive 
	//    z-direction with its up-vector in the y direction. The camera is translated
	//    by T and then rotated about the world origin by R. The C-arm system defines
	// 	  the rotation matricies using Euler angles with the convention 
	//
	//         R(α, β, γ) = Rz(α)Rx(β)Ry(γ), 
	//
	//    where Ri is a 3 × 3 matrix denoting rotation about the i-axis for i ∈{x, y, z}. 
	//    Here, α refers to the left-to-right anterior oblique rotational axis(LAO/RAO)
	//    and β refers to the cranial-caudal rotational axis(CRA/CAU), two common angles 
	//    in interventional radiology.
	//
	//    This yields a worldToCamera transform of |R~ | -t|, where R~ is the transpose 
	//    of cameraRotation and -t is the negative of the cameraTranslation vector and 
	//    |M | t| is a 3x4 matrix composed of 3x3 matrix M and 3x1 vector t. We 
	// 	  represent this matrix as a 4x4 matrix. 
	//
	//	  The 4x4 projection matrix from camera to image coordinates is 
	// 
	//							 |f 0 0 ox|, where f=focalLength and (ox, oy) is point 
	//    Mat4x4 cameraToImage = |0 f 0 oy|  on the detector closest to the camera. We
	//							 |0 0 f 0 |  assume this point is the center of the
	//							 |0 0 1 0 |  detector, i.e., (0, 0), in this implementation.
	// 
	//    The point (xv, yv, zv) in voxel coordinates is thus mapped along a camera ray 
	//    through (xv, yv, zv) to the homogeneous coordinate
	//    (xh, yh, zh, w) = cameraToImage * worldToCamera * voxelToWorld * (xv, yv, zv).
	// 
	// 	  The point on the detector (xd, yd, zd) is determined by applying the perspetive 
	//    divide to the homogenous image coordinate, i.e., (xd, yd, zd) = (xh/w, yh/w, zh/w).
	// 	  It can be shown that zd = f, i.e., the distance from the camera to the closest 
	//    point in the detector is the focal length f, as expected.
	//
	// 	  After applying the projection matrix, the detector point is mapped to pixel
	// 	  coordinates via the imageToPixel transform
	//
	//							|1/px 0    0    W/2|
	//	  Mat3x3 imageToPixel = |0    1/py 0    H/2|, where (px,py)=pixelSize, (W,H)=imageSize 
	//							|0    0    1    0  |
	//							|0    0    0    1  |
	//
	// 
	//    The pixel (i, j) is backprojected along the ray to the camera along the line from 
	//    (x, y, z) to the camera at t, where (x, y, z) = inv(imageToPixel) * (i, j, 0)
};