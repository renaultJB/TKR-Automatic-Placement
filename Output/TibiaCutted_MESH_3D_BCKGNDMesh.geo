Merge "Tibia.step";
Include "BCKGNDMesh.pos";

Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) splitmetis
Mesh.CharacteristicLengthFromPoints=0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthMin=1.0;
Mesh.CharacteristicLengthMax=4;
Mesh.Optimize = 1;

Background Mesh View[0];

SetOrder 2;
Mesh 1;
SetOrder 2;
Mesh 2;
SetOrder 2;
Mesh 3;
OptimizeMesh "Gmsh";
OptimizeMesh "Netgen";

SetOrder 2;

Save "Tibia.inp";

Exit;

