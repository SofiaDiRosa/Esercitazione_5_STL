#pragma once
#include <iostream>
#include "PolygonalMesh.hpp"

using namespace std;

namespace PolygonalLibrary
{
	// Importazione della PolygonalMesh e se ne testa la correttezza
	// mesh: PolygonalMesh struct
	// viene restituito il risultato della lettura: true se ha avuto successo, false altrimenti
	bool ImportMesh(PolygonalMesh& mesh);
	
	// Importazione delle proprietà di Cell0D dal file Cell0Ds.csv
	// mesh: PolygonalMesh struct
	// viene restituito il risultato della lettura: true se ha avuto successo, false altrimenti
	bool ImportCell0Ds(PolygonalMesh& mesh);
	
	// Importazione delle proprietà di Cell1D dal file Cell1Ds.csv
	// mesh: PolygonalMesh struct
	// viene restituito il risultato della lettura: true se ha avuto successo, false altrimenti
	bool ImportCell1Ds(PolygonalMesh& mesh);
	
	// Importazione delle proprietà di Cell2D dal file Cell2Ds.csv
	// mesh: PolygonalMesh struct
	// viene restituito il risultato della lettura: true se ha avuto successo, false altrimenti
	bool ImportCell2Ds(PolygonalMesh& mesh);
}