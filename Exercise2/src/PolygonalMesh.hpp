#pragma once
#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary {
	
struct PolygonalMesh
{
	unsigned int NumCell0Ds;
	vector<unsigned int> Cell0DsId = {}; //Id, dimensione: 1 x NumCell0Ds
	map<unsigned int, list<unsigned int>> Cell0DsMarker = {}; //Marker
	MatrixXd Cell0DsCoordinates = {}; //(X,Y), dimensione: 2 x NumCell0Ds
	
	unsigned int NumCell1Ds;
	vector<unsigned int> Cell1DsId = {}; //Id, dimensione: 1 x NumCell1Ds
	map<unsigned int, list<unsigned int>> Cell1DsMarker = {}; //Marker
	MatrixXi Cell1DsExtrems = {}; //(Origin, End), dimensione: 2 x NumCell1Ds
	
	unsigned int NumCell2Ds;
	vector<unsigned int> Cell2DsId = {}; //Id, dimensione: 1 x NumCell2Ds
	map<unsigned int, list<unsigned int>> Cell2DsMarker = {}; //Marker, ma è un vettore di 0
	vector<unsigned int> Cell2DsNumVertices = {}; //numero di vertici, dimensione: 1 x NumCell2Ds
	vector<array<unsigned int, 8>> Cell2DsVertices = {}; //vettore dei vertici di dimensione massima
	vector<unsigned int> Cell2DsNumEdges = {}; //numero di lati, dimensione: 1 x NumCell2Ds
	vector<array<unsigned int, 8>> Cell2DsEdges = {}; //vettore dei lati di dimensione massima
};

}
