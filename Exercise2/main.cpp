#include <iostream>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

bool TestEdges(const PolygonalMesh& mesh)
{
	const double epsilon = 1e-10;
	
	for(unsigned int cell_id = 0; cell_id < mesh.NumCell1Ds; cell_id++)
	{
		unsigned int origin_id = mesh.Cell1DsExtrems(0, cell_id);
		unsigned int end_id = mesh.Cell1DsExtrems(1, cell_id);
		
		Vector3d origin = mesh.Cell0DsCoordinates.col(origin_id);
		Vector3d end = mesh.Cell0DsCoordinates.col(end_id);
		
		double length = (end - origin).norm();
		if(length < epsilon)
		{
			cerr << "ERRORE: il segmento " << cell_id << " ha lunghezza nulla" << endl;
			return false;
		}
	}
	return true;
}

bool TestArea(const PolygonalMesh& mesh)
{
	const double epsilon = 1e-10;
	
	for(unsigned int cell_id = 0; cell_id < mesh.NumCell2Ds; cell_id ++)
	{
		const auto& vertices = mesh.Cell2DsVertices[cell_id];
		if(vertices.size() < 3)
			continue;
		
		// calcolo dell'area
		double area = 0.0;
		int n = vertices.size();
		for(int i = 0; i < n; i++)
		{
			unsigned int k = (i+1) % n;
			const auto& v1 = mesh.Cell0DsCoordinates.col(vertices[i]);
			const auto& v2 = mesh.Cell0DsCoordinates.col(vertices[k]);
			area += v1.x() * v2.y() - v1.y() * v2.x();
		}
		return abs(area) / 2.0;
		
		if(area < epsilon)
		{
			cerr << "ERRORE: il poligono con indice " << cell_id << " ha area nulla" << endl;
			return false;
		}
	}
	return true;
}

int main()
{
	PolygonalMesh mesh;
	
	if(!ImportMesh(mesh))
	{
		cerr << "File non trovato" << endl;
		return 1;
	}
	
	Gedim::UCDUtilities utilities;
	{
		vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);
		
		cell0Ds_properties[0].Label = "Marker";
		cell0Ds_properties[0].UnitLabel = "-";
		cell0Ds_properties[0].NumComponents = 1;
		
		vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
		for(const auto &m : mesh.Cell0DsMarker)
			for(const unsigned int id : m.second)
				cell0Ds_marker.at(id) = m.first;
			
		cell0Ds_properties[0].Data = cell0Ds_marker.data();
		utilities.ExportPoints("./Cell0Ds.inp",
							   mesh.Cell0DsCoordinates,
							   cell0Ds_properties);
	}
	
	{
		vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);
		
		cell1Ds_properties[0].Label = "Marker";
		cell1Ds_properties[0].UnitLabel = "-";
		cell1Ds_properties[0].NumComponents = 1;
		
		vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
		for(const auto &m : mesh.Cell1DsMarker)
			for(const unsigned int id : m.second)
				cell1Ds_marker.at(id) = m.first;
			
		cell1Ds_properties[0].Data = cell1Ds_marker.data();
		utilities.ExportSegments("./Cell1Ds.inp",
							   mesh.Cell0DsCoordinates,
							   mesh.Cell1DsExtrems,
							   {},
							   cell1Ds_properties);
	}
	
    return 0;
}
