Esize2D = ;

Merge "Tibia.stp";

Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) splitmetis
Mesh.CharacteristicLengthFromPoints=1;
Mesh.CharacteristicLengthMin=Esize2D;
Mesh.CharacteristicLengthMax=Esize2D;

Mesh 2;

Save "Tibia.msh";

Exit;