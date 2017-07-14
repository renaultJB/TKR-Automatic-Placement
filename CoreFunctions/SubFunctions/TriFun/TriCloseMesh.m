function [ TRout ] = TriCloseMesh( TRsup , TRin , nbElmts )
%TriCloseMesh : Perform a close morphological operation on a triangulation object 
% First a dilatation by nbElmts of TRin
% Second an erose by nbElmts of the newly created triangulation object TR

[ TR ] = TriDilateMesh( TRsup, TRin, nbElmts );
[ TRout ] = TriErodeMesh( TR, nbElmts );

end
