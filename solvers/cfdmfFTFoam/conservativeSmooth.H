// conservativeSmooth.H
#ifndef conservativeSmooth_H
#define conservativeSmooth_H

#include "fvCFD.H"

volVectorField conservativeSmooth
(
    const volVectorField& f,
    const fvMesh& mesh,
	const label nCFilterLoops
)
{
	volVectorField filteredField = f;
	for(int i = 1; i <= nCFilterLoops ; i++) //EA d2
	{
		// Smooth the field
		filteredField = fvc::average(fvc::interpolate(filteredField));
		filteredField.correctBoundaryConditions();
	}

	if (nCFilterLoops > 0)  // Skip if no filtering
	{
		// Compute conservation error
		vector netF_original = gSum(f.internalField() * mesh.V().field());
		vector netF_smoothed = gSum(filteredField.internalField() * mesh.V().field());
		vector error = netF_original - netF_smoothed;
		filteredField += dimensionedVector("fcorr", f.dimensions(), error / gSum(mesh.V()));
	}

    // Correct and return
    return filteredField;
}

#endif
