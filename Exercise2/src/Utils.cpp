#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace PolygonalLibrary
{
bool ImportMesh(PolygonalMesh& mesh)
{
	if(!ImportCell0Ds(mesh))
		return false;

	if(!ImportCell1Ds(mesh))
		return false;

	if(!ImportCell2Ds(mesh))
		return false;

	return true;
}
// ******************************************************************************
bool ImportCell0Ds(PolygonalMesh& mesh)
{
	ifstream file("./Cell0Ds.csv");
	
	if(file.fail())
		return false;
	
	list<vector<string>> listLines; //lista di vettori di stringhe
	
	string line;

	while(getline(file, line))
	{
		stringstream ss(line);
		string tmp;
		vector<string> val;
		while(getline(ss, tmp, ';'))
		{
			val.push_back(tmp);
		}
		listLines.push_back(val);
	}
	
	file.close();
	
	// remove header
	listLines.pop_front();
	
	mesh.NumCell0Ds = listLines.size();
	
	if(mesh.NumCell0Ds == 0)
	{
		cerr << "Non ci sono celle 0D" << endl;
		return false;
	}

	mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
	mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(2, mesh.NumCell0Ds);
	
	for(const vector<string>& line : listLines)
	{
		istringstream converter_id(line[0]);
		istringstream converter_marker(line[1]);
		istringstream converter_coordX(line[2]);
		istringstream converter_coordY(line[3]);
		
		unsigned int id;
		unsigned int marker;
		Vector2d coord;
		
		converter_id >> id;
		converter_marker >> marker;
		converter_coordX >> mesh.Cell0DsCoordinates(0, id);
		converter_coordY >> mesh.Cell0DsCoordinates(1, id);
		
		mesh.Cell0DsId.push_back(id);
		
		// memorizzazione dei marker
		if(marker != 0)
		{
			const auto it = mesh.Cell0DsMarker.find(marker);
			if(it == mesh.Cell0DsMarker.end())
			{
				mesh.Cell0DsMarker.insert({marker, {id}});
			}
			else
			{
				mesh.Cell0DsMarker[marker].push_back(id);
			}
		}
	}
	
	return true;
}
// ******************************************************************************
bool ImportCell1Ds(PolygonalMesh& mesh)
{
	ifstream file("./Cell1Ds.csv");
	
	if(file.fail())
		return false;
	
	list<vector<string>> listLines; //lista di vettori di stringhe
	
	string line;

	while(getline(file, line))
	{
		stringstream ss(line);
		string tmp;
		vector<string> val;
		while(getline(ss, tmp, ';'))
		{
			val.push_back(tmp);
		}
		listLines.push_back(val);
	}
	
	file.close();
	
	// remove header
	listLines.pop_front();
	
	mesh.NumCell1Ds = listLines.size();
	
	if(mesh.NumCell1Ds == 0)
	{
		cerr << "Non ci sono celle 1D" << endl;
		return false;
	}
	
	mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
	mesh.Cell1DsExtrems = Eigen::MatrixXi(2, mesh.NumCell1Ds);
	
	for(const vector<string>& line : listLines)
	{
		istringstream converter_id(line[0]);
		istringstream converter_marker(line[1]);
		istringstream converter_origin(line[2]);
		istringstream converter_end(line[3]);
		unsigned int id;
		unsigned int marker;
		Vector2i vertices;
		
		converter_id >> id;
		converter_marker >> marker;
		converter_origin >> mesh.Cell1DsExtrems(0, id);
		converter_end >> mesh.Cell1DsExtrems(1, id);
		
		mesh.Cell1DsId.push_back(id);
		
		// memorizzazione dei marker
		if(marker != 0)
		{
			const auto it = mesh.Cell1DsMarker.find(marker);
			if(it == mesh.Cell1DsMarker.end())
			{
				mesh.Cell1DsMarker.insert({marker, {id}});
			}
			else
			{
				mesh.Cell1DsMarker[marker].push_back(id);
			}
		}
	}
	
	return true;
}
// ******************************************************************************
bool ImportCell2Ds(PolygonalMesh& mesh)
{
	ifstream file("./Cell2Ds.csv");
	
	if(file.fail())
		return false;
	
	list<vector<string>> listLines; //lista di vettori di stringhe
	
	string line;

	while(getline(file, line))
	{
		stringstream ss(line);
		string tmp;
		vector<string> val;
		while(getline(ss, tmp, ';'))
		{
			val.push_back(tmp);
		}
		listLines.push_back(val);
	}
	
	file.close();
	
	// remove header
	listLines.pop_front();
	
	mesh.NumCell2Ds = listLines.size();
	
	if(mesh.NumCell2Ds == 0)
	{
		cerr << "Non ci sono celle 2D" << endl;
		return false;
	}
	
	mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
	mesh.Cell2DsNumVertices.reserve(mesh.NumCell2Ds);
	mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
	mesh.Cell2DsNumEdges.reserve(mesh.NumCell2Ds);
	mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);
	
	for(const vector<string>& line : listLines)
	{
		istringstream converter_id(line[0]);
		istringstream converter_marker(line[1]);
		istringstream converter_numver(line[2]);
		istringstream converter_vertices(line[3]);
		istringstream converter_numedg(line[4]);
		istringstream converter_edges(line[5]);
		unsigned int id;
		unsigned int marker;
		unsigned int numver;
		array<unsigned int, 8> vertices;
		unsigned int numedg;
		array<unsigned int, 8> edges;
		
		converter_id >> id;
		converter_marker >> marker;
		converter_numver >> numver;
		for(unsigned int i = 3; i < 8; i++)
			converter_vertices >> vertices[i];
		converter_numedg >> numedg;
		for(unsigned int i = 3; i < 8; i++)
			converter_edges >> edges[i];
		
		mesh.Cell2DsId.push_back(id);
		mesh.Cell2DsNumVertices.push_back(numver);
		mesh.Cell2DsVertices.push_back(vertices);
		mesh.Cell2DsNumEdges.push_back(numedg);
		mesh.Cell2DsEdges.push_back(edges);
		
		// memorizzazione dei marker
		if(marker != 0)
		{
			const auto it = mesh.Cell2DsMarker.find(marker);
			if(it == mesh.Cell2DsMarker.end())
			{
				mesh.Cell2DsMarker.insert({marker, {id}});
			}
			else
			{
				mesh.Cell2DsMarker[marker].push_back(id);
			}
		}
	}
	
	return true;
}

}