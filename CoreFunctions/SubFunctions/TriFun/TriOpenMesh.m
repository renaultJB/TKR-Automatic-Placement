function [ TRout ] = TriOpenMesh( TRsup , TRin , nbElmts )
%TriOpenMesh : Perform a "open" morphological operation on a triangulation object 
% First a erosion by nbElmts of TRin
% Second an dilatation by nbElmts of the newly created triangulation object TR


[ TR ] = TriErodeMesh( TRin, nbElmts );
[ TRout ] = TriDilateMesh( TRsup, TR, nbElmts );

end

