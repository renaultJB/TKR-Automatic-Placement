function [ TRout ] = TriKeepLargestPatch( TRin )
%TriKeepLargestPatch : Keep the largest (by area) connected patch of a triangulation
% Object

% Unconnect patch potentially sharing only one node between 
Trin2 = TriErodeMesh(TRin,1);

% Get the exterior list of edges on the boundary of the triangulation
Segments = Trin2.freeBoundary;

% Recreate the the orientated close(s) border curves
% 
j=1;
Curves=struct();
i=1;
while ~isempty(Segments)
    Curves(i).NodesID = zeros(length(Segments),1);
    Curves(i).NodesID(j)=Segments(1,1);
    Curves(i).NodesID(j+1)=Segments(1,2);
    Segments(1,:)=[];
    j=j+1;
    [Is,Js] = ind2sub(size(Segments),find(Segments(:) == Curves(i).NodesID(j)));
    Nk = Segments(Is,round(Js+2*(1.5-Js)));
    Segments(Is,:)=[];
    j=j+1;
    while ~isempty(Nk)
        Curves(i).NodesID(j) = Nk(1);
        [Is,Js] = ind2sub(size(Segments),find(Segments(:) == Curves(i).NodesID(j)));
        Nk = Segments(Is,round(Js+2*(1.5-Js)));
        Segments(Is,:)=[];
        j=j+1;
    end
    Curves(i).NodesID(Curves(i).NodesID==0) = []  ;
    Curves(i).NodesID;
    i=i+1;
end


% Identify the largest patch of the trinagulation
if length(Curves)>1
    
    SizePatch = zeros(length(Curves),1);
    Patchs = struct();
    for k = 1 : length(Curves)
        Patchs(k).TR = TriConnectedPatch(TRin,Trin2.Points(Curves(k).NodesID(:),:));
        [ Properties ] = TriMesh2DProperties( Patchs(k).TR );
        SizePatch(k) = Properties.TotalArea;
    end
    [~,IMax] = max(SizePatch);
    TRout = Patchs(IMax).TR ;
    
else
    TRout = TriDilateMesh(TRin,Trin2,1);
end



end

