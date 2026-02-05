/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  9.0                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     9.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)

//*******************************************************************************
m4_define(xLength, 0.031) // the width of straight Channel
m4_define(yLength, 0.031) // the width of straight Channel
m4_define(zLength, 0.350) // the height of straight Channel

//mesh
m4_define(NmeshX, 31) //number of radial mesh
m4_define(NmeshY, 31) //number of radial mesh
m4_define(NmeshZ, 270) //number of radial mesh

convertToMeters 1;
vertices
(
//**********************************buttom nodes*********************************
(-calc(xLength/2) -calc(yLength/2) 0) vlabel(B1)
(calc(xLength/2) -calc(yLength/2) 0) vlabel(B2)
(calc(xLength/2) calc(yLength/2) 0) vlabel(B3)
(-calc(xLength/2) calc(yLength/2) 0) vlabel(B4)

//********************************** upper nodes *********************************
(-calc(xLength/2) -calc(yLength/2) zLength) vlabel(U1)
(calc(xLength/2) -calc(yLength/2) zLength) vlabel(U2)
(calc(xLength/2) calc(yLength/2) zLength) vlabel(U3)
(-calc(xLength/2) calc(yLength/2) zLength) vlabel(U4)

//************************************end upper nodes*****************************
);

// Defining blocks:

blocks
(

    // block1 
    hex (B1 B2 B3 B4 U1 U2 U3 U4) AB
        (NmeshX NmeshY NmeshZ)
    simpleGrading (1 1 1)
);

edges
(
  //line B1 U1 
  //line B2 U2 
  //line B3 U3 
  //line B4 U4

  //line B1 B2 
  //line B2 B3 
  //line B3 B4 
  //line B4 B1

  //line U1 U2 
  //line U2 U3 
  //line U3 U4 
  //line U4 U1  
);

boundary
(

    inlet
    {
		type patch;
		faces
		(
		   (B1 B2 B3 B4)
		);
    }

    walls
    {
		type patch;
		faces
		(
		   (B1 B2 U2 U1)
		   (B2 B3 U3 U2)
		   (B3 B4 U4 U3)
		   (B4 B1 U1 U4)
		);  
    }

    outlet
    {
		type patch;
		faces
		(
		   (U1 U2 U3 U4)	   
		);
    }
);

mergePatchPairs 
(
);
