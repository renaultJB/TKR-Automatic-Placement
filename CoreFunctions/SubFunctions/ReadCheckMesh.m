function [ TR1, TR2 ] = ReadCheckMesh( MeshFileName1, MeshFileName2 )
%ReadAndCheckMESH Read data from GMSH 2D Mesh Files and Parse to matlab then
%check that the normal are correctly oriented (Pointing Outward)
% INPUT(s)
% MeshFileName1 : Char list  of Fullname -> Path + name of the file 1
% Optional 2nd Read mesh :
% MeshFileName2 : Char list  of Fullname -> Path + name of the file 2
% OUTPUT(s) :
% TR1 : Triangulation Object of the file 1
% Optional if 2nd Read mesh :
% TR1 : Triangulation Object  of the file 1



XYZELMTS = py.txt2mtlb.read_meshGMSH(MeshFileName1);
Pts2D = [cell2mat(cell(XYZELMTS{'X'}))' cell2mat(cell(XYZELMTS{'Y'}))' cell2mat(cell(XYZELMTS{'Z'}))'];
Elmts2D = double([cell2mat(cell(XYZELMTS{'N1'}))' cell2mat(cell(XYZELMTS{'N2'}))' cell2mat(cell(XYZELMTS{'N3'}))']);


% Verify that normal are outward-pointing and fix if not
Elmts2D = fixNormals( Pts2D, Elmts2D );
TR1 = triangulation(Elmts2D,Pts2D);

if nargin >1 && nargout==nargin
    %Read distal Tibia
    XYZELMTS = py.txt2mtlb.read_meshGMSH(MeshFileName2);
    Pts2D = [cell2mat(cell(XYZELMTS{'X'}))' cell2mat(cell(XYZELMTS{'Y'}))' cell2mat(cell(XYZELMTS{'Z'}))'];
    Elmts2D = double([cell2mat(cell(XYZELMTS{'N1'}))' cell2mat(cell(XYZELMTS{'N2'}))' cell2mat(cell(XYZELMTS{'N3'}))']);

    % Verify that normal are outward-pointing and fix if not
    Elmts2D = fixNormals( Pts2D, Elmts2D );
    TR2 = triangulation(Elmts2D,Pts2D);
elseif nargin+nargout > 2
    msg = 'Error : Number of inputs and outputs are not consistent';
    error(msg)   
end


end

