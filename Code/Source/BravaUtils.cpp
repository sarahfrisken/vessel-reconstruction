//
// BravaUtils.cpp
// Implementation of BravaUtils.
//

#include "BravaUtils.h"

#include <QVector3D>

// 
// Public
//
BravaUtils::BravaUtils()
{

}
BravaUtils::~BravaUtils()
{

}

// 
// Public
//
void BravaUtils::importBravaData(const char* filename)
{
	// Open the file
	FILE* fp;
	fopen_s(&fp, filename, "r");
	if (!fp) return;

	// Read the file line by line
	int lineRead = 1;
	std::vector<VesselModel::Vertex> inputVertices;

	do {
		char line[100];
		lineRead = (fgets(line, 100, fp) == NULL) ? 0 : 1;
		if (lineRead && line[0] != '#') {
			int ID;
			int branchIdx;
			int parentIdx;
			QVector3D pos;
			float radius;
			float time;
			QVector3D velocity;
			int numRead = sscanf_s(line, "%d %d %f %f %f %f %d",
				&ID, &branchIdx, &(pos[0]), &(pos[2]), &(pos[1]), &radius, &parentIdx);
			if (numRead != 7) exit(0);

			if (m_branches.size() == 0) {
				// Create the first branch
				Branch branch;
				branch.ID = branchIdx;
				m_branches.push_back(branch);
				time = 0;
				velocity = {0, 0, 0}; // Not ideal but less complicated than finding next vertex
			}
			else if ((ID != parentIdx + 1) || (branchIdx != m_branches.back().ID)) {
				// Begin a new branch. Duplicate the parent node and add it to the new branch.
				VesselModel::Vertex first = inputVertices[size_t(parentIdx) - 1];
				Branch branch;
				branch.ID = branchIdx;
				branch.vertices.push_back(first);
				m_branches.push_back(branch);
				QVector3D vecToPos;
				vecToPos = pos - first.pos;
				time = first.time + vecToPos.length();
				velocity = first.velocity;
			}
			else {
				VesselModel::Vertex prev = m_branches.back().vertices.back();
				QVector3D vecFromPrev;
				vecFromPrev = pos - prev.pos;
				time = m_branches.back().vertices.back().time + vecFromPrev.length(); // Assume constant normalized speed: time propotional to distance from root
				velocity = vecFromPrev.normalized(); // Assume constant speed: velocity is in direction of branch
			}

			// Add the new point to the current branch and the temporary list of input points
			VesselModel::Vertex vertex = { pos, radius, time, velocity }; 
			inputVertices.push_back(vertex);
			m_branches.back().vertices.push_back(vertex);
		}
	} while (lineRead);

	// Close the file
	fclose(fp);
}

//
// Private
//
std::vector<VesselModel::Line> BravaUtils::vesselLines(std::vector<int> branchIDs)
{
	// Returns vessel segments of the specified branches. Each vessel segment is a line
	// that has two vertices (endpoints).
	std::vector<VesselModel::Line> lines;
	for (std::vector<int>::iterator id = branchIDs.begin(); id != branchIDs.end(); ++id) {
		for (std::vector<Branch>::iterator b = m_branches.begin(); b != m_branches.end(); ++b) {
			if (b->ID != *id) continue;
			VesselModel::Vertex prevV = b->vertices[0];
			for (int i = 1; i < b->vertices.size(); i++) {
				VesselModel::Vertex v = b->vertices[i];
				lines.push_back(VesselModel::Line({prevV, v}));
				prevV = v;
			}
		}
	}
	return lines;
}

