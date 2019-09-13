Merge "Tibia.step";
Mesh.RemeshAlgorithm=1; //(0) nosplit (1) automatic (2) splitmetis
Mesh.CharacteristicLengthFromPoints=1;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.CharacteristicLengthMin=1.25;
Mesh.CharacteristicLengthMax=4;
Mesh.Optimize = 1;
Mesh.OptimizeNetGen = 1;
Mesh.HighOrderOptimize = 1;


Zcent = ZALT; // Altitude de la bande raffinée
Largeur = 15 ; // Demi Largeur de la bande raffinée
MinElSize = 2.5 ;
MaxElSize = 3.5 ;
d = 3.2; // Coefficient de largeur de la transition au bord de la boxe entre 0 et 10

a = (MaxElSize-MinElSize)/((d^2-1)*Largeur^2);
b = (d^2*MinElSize-MaxElSize)/(d^2-1);

Field[1] = Box;
Field[1].VIn = MinElSize;
Field[1].VOut = MaxElSize;
Field[1].XMin = -1000; 
Field[1].XMax = 1000;
Field[1].YMin = -1000;
Field[1].YMax = 1000;
Field[1].ZMin = Zcent - Largeur;
Field[1].ZMax = Zcent + Largeur;

Field[2] = MathEval;
Field[2].F = Sprintf("Fabs(%+g*(z%+g)^2%+g%+g)%+g",a,-Zcent,b,-MinElSize,MinElSize);
Printf("Fabs(%+g*(z%+g)^2%+g%+g)%+g",a,-Zcent,b,-MinElSize,MinElSize);

Field[3] = Min;
Field[3].FieldsList = {1,2};
Background Field = 3;

SetOrder 2;

Mesh 2;
Mesh 3;

OptimizeMesh "Gmsh";
OptimizeMesh "Netgen";

SetOrder 2;

Mesh.HighOrderOptimize = 1;


Save "Tibia.inp";

Exit;

