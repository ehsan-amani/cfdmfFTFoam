clear
clc
%% This is a script to convert BINARY STL mesh (from Salome or any other Program) to
%% a mesh format readable by cfdmfFTFoam solver as an input mesh. 


%%-------Change These Based on phase properties
    bubbleNo = 1;                   %No of bubbles in domain
    outerFluidDensity = 1240;       %Density of Continuous  phase
    density = 955;                  %Density of bubble/drop
    outerFluidViscosity = 262e-03;  %Viscosity of Continuous  phase
    viscosity = 29e-03;             %Viscosity of bubble/drop
    surfaceTensionCoeff = 29e-03;   %Surface Tension Coefficient
    Eotvos = 0;
    Morton = 0;
    diameter = 0.038404;            %Sphere Equivalent diameter of bubble
    convertToMeters = 0.001;        %Convert Salome mesh from mm to meters
%%--------------------------------------------------



%%--------Specify the name of your STL file here: 
inputName = 'taylorDrop.stl';

%%--------Specify the output filename here:
outputName = 'initialFrontMesh';

%% After running this code a file named "initialFrontMesh" would be created
%% Copy this file to FTResult directory of your case.
[v, f, n, c, stltitle] = stlread(inputName);

convFTmsh(outputName,f,v,bubbleNo, outerFluidDensity, ...
   density,  outerFluidViscosity, viscosity, surfaceTensionCoeff, Eotvos, ...
   Morton,diameter,convertToMeters);
