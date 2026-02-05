/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// 090512: Original files by Håkan Nilsson, Chalmers University of Technology
// 2014   Modified by Ehsan Amani, Amirkabir University of Technology.

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

//Geometry
m4_define(rIn, 1) //Small radius ok
//m4_define(rOut, 100) //Large radius ok
//m4_define(chLength, 500) //Chamber lenght ok
m4_define(cin, 6) //inlet lenght ok

//Number of cells:
//Radial direction:
m4_define(rNumberOfCells1st, 15)  // in the first O-grid
//m4_define(rNumberOfCells2nd, 30)  // in the second O-grid
//m4_define(rNumberOfCells3rd, 10)  // in the third O-grid
//Tangential direction:
m4_define(tNumberOfCells, 38)
//Axial direction:
m4_define(zABnumberOfCells, 200)
//m4_define(zBCnumberOfCells, 85)

//Radial grading:
m4_define(rGrading1, 1)
//m4_define(rGrading2, 8)
//m4_define(rGrading3, 1)

//Axial grading:
m4_define(zGradingAB, 1)
//m4_define(zGradingBC, 8)

//*******************************************************************************

//Plane A:
m4_define(zA, 0)
m4_define(rA, rIn)
m4_define(rRelA, 0.7)
m4_define(rRelAc, 0.8)

//Plane B:
m4_define(zB, cin)
//m4_define(rB, rOut)
m4_define(rB1st, rIn)
//m4_define(rB2nd, calc(0.5*(rB1st+rB)))
m4_define(rRelB, 0.7)
m4_define(rRelBc, 0.8)

//Plane C:
//m4_define(zC, chLength)
//m4_define(rC, rOut)
//m4_define(rC1st, rIn)
//m4_define(rC2nd, calc(0.5*(rC1st+rC)))
//m4_define(rRelC, 0.7)
//m4_define(rRelCc, 0.8)

//*******************************************************************************

convertToMeters 1;

vertices
(
//Plane A:
(calc(rRelA*rA*cos(pi/4)) -calc(rRelA*rA*sin(pi/4)) zA) vlabel(A0)
(calc(rRelA*rA*cos(pi/4)) calc(rRelA*rA*sin(pi/4)) zA) vlabel(A1)
(calc(-rRelA*rA*cos(pi/4)) calc(rRelA*rA*sin(pi/4)) zA) vlabel(A2)
(calc(-rRelA*rA*cos(pi/4)) -calc(rRelA*rA*sin(pi/4)) zA) vlabel(A3)
(calc(rA*cos(pi/4)) -calc(rA*sin(pi/4)) zA) vlabel(A4)
(calc(rA*cos(pi/4)) calc(rA*sin(pi/4)) zA) vlabel(A5)
(calc(-rA*cos(pi/4)) calc(rA*sin(pi/4)) zA) vlabel(A6)
(calc(-rA*cos(pi/4)) -calc(rA*sin(pi/4)) zA) vlabel(A7)

//Plane B:
(calc(rRelB*rB1st*cos(pi/4)) -calc(rRelB*rB1st*sin(pi/4)) zB) vlabel(B0)
(calc(rRelB*rB1st*cos(pi/4)) calc(rRelB*rB1st*sin(pi/4)) zB) vlabel(B1)
(calc(-rRelB*rB1st*cos(pi/4)) calc(rRelB*rB1st*sin(pi/4)) zB) vlabel(B2)
(calc(-rRelB*rB1st*cos(pi/4)) -calc(rRelB*rB1st*sin(pi/4)) zB) vlabel(B3)
(calc(rB1st*cos(pi/4)) -calc(rB1st*sin(pi/4)) zB) vlabel(B4)
(calc(rB1st*cos(pi/4)) calc(rB1st*sin(pi/4)) zB) vlabel(B5)
(calc(-rB1st*cos(pi/4)) calc(rB1st*sin(pi/4)) zB) vlabel(B6)
(calc(-rB1st*cos(pi/4)) -calc(rB1st*sin(pi/4)) zB) vlabel(B7)
//(calc(rB2nd*cos(pi/4)) -calc(rB2nd*sin(pi/4)) zB) vlabel(B8)
//(calc(rB2nd*cos(pi/4)) calc(rB2nd*sin(pi/4)) zB) vlabel(B9)
//(calc(-rB2nd*cos(pi/4)) calc(rB2nd*sin(pi/4)) zB) vlabel(B10)
//(calc(-rB2nd*cos(pi/4)) -calc(rB2nd*sin(pi/4)) zB) vlabel(B11)
//(calc(rB*cos(pi/4)) -calc(rB*sin(pi/4)) zB) vlabel(B8)
//(calc(rB*cos(pi/4)) calc(rB*sin(pi/4)) zB) vlabel(B9)
//(calc(-rB*cos(pi/4)) calc(rB*sin(pi/4)) zB) vlabel(B10)
//(calc(-rB*cos(pi/4)) -calc(rB*sin(pi/4)) zB) vlabel(B11)
//(calc(rB*cos(pi/4)) -calc(rB*sin(pi/4)) zB) vlabel(B12)
//(calc(rB*cos(pi/4)) calc(rB*sin(pi/4)) zB) vlabel(B13)
//(calc(-rB*cos(pi/4)) calc(rB*sin(pi/4)) zB) vlabel(B14)
//(calc(-rB*cos(pi/4)) -calc(rB*sin(pi/4)) zB) vlabel(B15)

//Plane C:
//(calc(rRelC*rC1st*cos(pi/4)) -calc(rRelC*rC1st*sin(pi/4)) zC) vlabel(C0)
//(calc(rRelC*rC1st*cos(pi/4)) calc(rRelC*rC1st*sin(pi/4)) zC) vlabel(C1)
//(calc(-rRelC*rC1st*cos(pi/4)) calc(rRelC*rC1st*sin(pi/4)) zC) vlabel(C2)
//(calc(-rRelC*rC1st*cos(pi/4)) -calc(rRelC*rC1st*sin(pi/4)) zC) vlabel(C3)
//(calc(rC1st*cos(pi/4)) -calc(rC1st*sin(pi/4)) zC) vlabel(C4)
//(calc(rC1st*cos(pi/4)) calc(rC1st*sin(pi/4)) zC) vlabel(C5)
//(calc(-rC1st*cos(pi/4)) calc(rC1st*sin(pi/4)) zC) vlabel(C6)
//(calc(-rC1st*cos(pi/4)) -calc(rC1st*sin(pi/4)) zC) vlabel(C7)
//(calc(rC2nd*cos(pi/4)) -calc(rC2nd*sin(pi/4)) zC) vlabel(C8)
//(calc(rC2nd*cos(pi/4)) calc(rC2nd*sin(pi/4)) zC) vlabel(C9)
//(calc(-rC2nd*cos(pi/4)) calc(rC2nd*sin(pi/4)) zC) vlabel(C10)
//(calc(-rC2nd*cos(pi/4)) -calc(rC2nd*sin(pi/4)) zC) vlabel(C11)
//(calc(rC*cos(pi/4)) -calc(rC*sin(pi/4)) zC) vlabel(C8)
//(calc(rC*cos(pi/4)) calc(rC*sin(pi/4)) zC) vlabel(C9)
//(calc(-rC*cos(pi/4)) calc(rC*sin(pi/4)) zC) vlabel(C10)
//(calc(-rC*cos(pi/4)) -calc(rC*sin(pi/4)) zC) vlabel(C11)
//(calc(rC*cos(pi/4)) -calc(rC*sin(pi/4)) zC) vlabel(C12)
//(calc(rC*cos(pi/4)) calc(rC*sin(pi/4)) zC) vlabel(C13)
//(calc(-rC*cos(pi/4)) calc(rC*sin(pi/4)) zC) vlabel(C14)
//(calc(-rC*cos(pi/4)) -calc(rC*sin(pi/4)) zC) vlabel(C15)
);


// Defining blocks:
blocks
(

    //Blocks between plane A and plane B:
    // block0 - positive x O-grid block
    hex (A5 A1 A0 A4 B5 B1 B0 B4 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 zGradingAB)
    // block1 - positive y O-grid block
    hex (A6 A2 A1 A5 B6 B2 B1 B5 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 zGradingAB)
    // block2 - negative x O-grid block
    hex (A7 A3 A2 A6 B7 B3 B2 B6 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 zGradingAB)
    // block3 - negative y O-grid block
    hex (A4 A0 A3 A7 B4 B0 B3 B7 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 zGradingAB)
    // block4 - central O-grid block
    hex (A0 A1 A2 A3 B0 B1 B2 B3 ) AB
    (tNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (1 1 zGradingAB)

/*
    //Blocks between plane B and plane C:
    // block0 - positive x O-grid block 1st belt
    hex (B5 B1 B0 B4 C5 C1 C0 C4 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 zGradingBC)
    // block1 - positive y O-grid block 1st belt
    hex (B6 B2 B1 B5 C6 C2 C1 C5 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 zGradingBC)
    // block2 - negative x O-grid block 1st belt
    hex (B7 B3 B2 B6 C7 C3 C2 C6 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 zGradingBC)
    // block3 - negative y O-grid block 1st belt
    hex (B4 B0 B3 B7 C4 C0 C3 C7 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 zGradingBC)
    // block4 - central O-grid block 
    hex (B0 B1 B2 B3 C0 C1 C2 C3 ) BC
    (tNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (1 1 zGradingBC)
    // block5 - positive x O-grid block 2nd belt 
    hex (B9 B5 B4 B8 C9 C5 C4 C8 ) BC
    (rNumberOfCells2nd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 zGradingBC)
    // block6 - positive y O-grid block 2nd belt
    hex (B10 B6 B5 B9 C10 C6 C5 C9 ) BC
    (rNumberOfCells2nd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 zGradingBC)
    // block7 - negative x O-grid block 2nd belt
    hex (B11 B7 B6 B10 C11 C7 C6 C10 ) BC
    (rNumberOfCells2nd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 zGradingBC)
    // block8 - negative y O-grid block 2nd belt
    hex (B8 B4 B7 B11 C8 C4 C7 C11 ) BC
    (rNumberOfCells2nd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 zGradingBC)

    // block9 - positive x O-grid block 3rd belt
    hex (B13 B9 B8 B12 C13 C9 C8 C12 ) BC
    (rNumberOfCells3rd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading3 1 zGradingBC)
    // block10 - positive y O-grid block 3rd belt
    hex (B14 B10 B9 B13 C14 C10 C9 C13 ) BC
    (rNumberOfCells3rd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading3 1 zGradingBC)
    // block11 - negative x O-grid block 3rd belt
    hex (B15 B11 B10 B14 C15 C11 C10 C14 ) BC
    (rNumberOfCells3rd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading3 1 zGradingBC)
    // block12 - negative y O-grid block 3rd belt
    hex (B12 B8 B11 B15 C12 C8 C11 C15 ) BC
    (rNumberOfCells3rd tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading3 1 zGradingBC)
*/
);

edges
(
    
    //Plane A:
    arc A0 A1 (calc(rRelAc*rRelA*rA) 0 zA)
    arc A1 A2 (0 calc(rRelAc*rRelA*rA) zA)
    arc A2 A3 (-calc(rRelAc*rRelA*rA) 0 zA)
    arc A3 A0 (0 -calc(rRelAc*rRelA*rA) zA)
    arc A4 A5 (rA 0 zA)
    arc A5 A6 (0 rA zA)
    arc A6 A7 (-rA 0 zA)
    arc A7 A4 (0 -rA zA)

    //Plane B:
    arc B0  B1  (calc(rRelBc*rRelB*rB1st) 0 zB)
    arc B1  B2  (0 calc(rRelBc*rRelB*rB1st) zB)
    arc B2  B3  (-calc(rRelBc*rRelB*rB1st) 0 zB)
    arc B3  B0  (0 -calc(rRelBc*rRelB*rB1st) zB)
    arc B4  B5  (rB1st 0 zB)
    arc B5  B6  (0 rB1st zB)
    arc B6  B7  (-rB1st 0 zB)
    arc B7  B4  (0 -rB1st zB)
    //arc B8  B9  (rB2nd 0 zB)
    //arc B9  B10 (0 rB2nd zB)
    //arc B10 B11 (-rB2nd 0 zB)
    //arc B11 B8  (0 -rB2nd zB)
    //arc B8  B9  (rB 0 zB)
    //arc B9  B10 (0 rB zB)
    //arc B10 B11 (-rB 0 zB)
    //arc B11 B8  (0 -rB zB)
    //arc B12 B13 (rB 0 zB)
    //arc B13 B14 (0 rB zB)
    //arc B14 B15 (-rB 0 zB)
    //arc B15 B12 (0 -rB zB)

    //Plane C:
    //arc C0  C1 (calc(rRelCc*rRelC*rC1st) 0 zC)
    //arc C1  C2 (0 calc(rRelCc*rRelC*rC1st) zC)
    //arc C2  C3 (-calc(rRelCc*rRelC*rC1st) 0 zC)
    //arc C3  C0 (0 -calc(rRelCc*rRelC*rC1st) zC)
    //arc C4  C5 (rC1st 0 zC)
    //arc C5  C6 (0 rC1st zC)
    //arc C6  C7 (-rC1st 0 zC)
    //arc C7  C4 (0 -rC1st zC)
    //arc C8  C9  (rC2nd 0 zC)
    //arc C9  C10 (0 rC2nd zC)
    //arc C10 C11 (-rC2nd 0 zC)
    //arc C11 C8  (0 -rC2nd zC)
    //arc C8  C9  (rC 0 zC)
    //arc C9  C10 (0 rC zC)
    //arc C10 C11 (-rC 0 zC)
    //arc C11 C8  (0 -rC zC)
    //arc C12 C13 (rC 0 zC)
    //arc C13 C14 (0 rC zC)
    //arc C14 C15 (-rC 0 zC)
    //arc C15 C12 (0 -rC zC)

);

// Defining patches:
boundary
(
    inlet
    {
		type patch;
		faces
		(
		   (A1 A5 A4 A0)
		   (A2 A6 A5 A1)
		   (A3 A7 A6 A2)
		   (A0 A4 A7 A3)
		   (A3 A2 A1 A0)
		);
    }
    walls
    {
		type wall;
		faces
		(
		   (A4 A5 B5 B4)
		   (A5 A6 B6 B5)
		   (A6 A7 B7 B6)
		   (A7 A4 B4 B7)

		   //(B5 B9 B8 B4)
		   //(B6 B10 B9 B5)
		   //(B7 B11 B10 B6)
		   //(B4 B8 B11 B7)
		   //(B9 B13 B12 B8)
		   //(B10 B14 B13 B9)
		   //(B11 B15 B14 B10)
		   //(B8 B12 B15 B11)

		   //(B12 B13 C13 C12)
		   //(B13 B14 C14 C13)
		   //(B14 B15 C15 C14)
		   //(B15 B12 C12 C15)
		   //(B8 B9 C9 C8)
		   //(B9 B10 C10 C9)
		   //(B10 B11 C11 C10)
		   //(B11 B8 C8 C11)   
		);  
    }
    outlet
    {
		type patch;
		faces
		(
		   //(C0 C4 C5 C1)
		   //(C1 C5 C6 C2)
		   //(C2 C6 C7 C3)
		   //(C3 C7 C4 C0)
		   //(C0 C1 C2 C3)
		   (B0 B4 B5 B1)
		   (B1 B5 B6 B2)
		   (B2 B6 B7 B3)
		   (B3 B7 B4 B0)
		   (B0 B1 B2 B3)

		   //(C4 C8 C9 C5)
		   //(C5 C9 C10 C6)
		   //(C6 C10 C11 C7)
		   //(C7 C11 C8 C4)

		   //(C8 C12 C13 C9)
		   //(C9 C13 C14 C10)
		   //(C10 C14 C15 C11)
		   //(C11 C15 C12 C8)
		);
    }
);

mergePatchPairs 
(
);

// ************************************************************************* //
