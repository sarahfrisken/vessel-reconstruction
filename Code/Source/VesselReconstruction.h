//
// VesselReconstruction.h
// User interface for this vessel reconstruction project.
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

#include <QtWidgets/QMainWindow>
#include "ui_VesselReconstruction.h"

class VesselReconstruction : public QMainWindow
{
    Q_OBJECT

public:
    VesselReconstruction(QWidget *parent = nullptr);
    ~VesselReconstruction();

private slots:
    void onGenerateXmlTemplates();
    void onPreprocessBravaModels();
    void onReconstruct();
    void onValidate();

private:
    Ui::VesselReconstructionClass ui;
    float m_voxelSize;
    float m_maxDist;
};
