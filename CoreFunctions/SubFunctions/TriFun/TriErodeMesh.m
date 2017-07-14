function [ TRout ] = TriErodeMesh( TRin, nbElmts )
% Erode a triangulation TRin :
%   remove elements from TRin triangulation that have an edge   
%   on the the border of TRin
% Update TRin and Repeat the operation nbElmts times

BorderNodesID = unique(TRin.freeBoundary);
ElmtsInitial = TRin.ConnectivityList;
ElmtsBorder = find(sum(ismember(ElmtsInitial,BorderNodesID),2)>0);
ElmtsInitial = ElmtsBorder;

if nbElmts>1
    for i = 1 : nbElmts-1
        ElmtNeighbours = unique(NotNaN(TRin.neighbors(ElmtsInitial)));
        ElmtsInitial = unique(ElmtNeighbours);
        ElmtsBorder = unique([ ElmtsBorder ; ElmtsInitial ]);
    end
end

ElmtsKept = ones(length(TRin.ConnectivityList),1);
ElmtsKept(ElmtsBorder) = 0;

TRout = TriReduceMesh( TRin, find(ElmtsKept));

end
