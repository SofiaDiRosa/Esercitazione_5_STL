#include <iostream>
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

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
