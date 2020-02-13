function [ TRout, HolePtsProj ] = TriFillPlanarHoles( TRin, moveNodes)

%TRIFILLPLANARHOLES Fill planar convex holes in the triangulation
%   For now the holes have to be planar
%   FOR NOW WORKS WITH ONLY ONE HOLE

FB = TRin.freeBoundary;
Nodes0 = TRin.Points;



% Get the least square plan of the hole
HolePts = Nodes0(FB(:,1),:);
[V,~] = eig(cov(HolePts)); n = V(:,1);
HoleCenter = mean(HolePts);
TriCenter =  mean(Nodes0);



if nargin>1 && moveNodes

    % Project the Edge vertices on the plan
    moveNodes
    HolePtsProj = HolePts - repmat(bsxfun(@minus,HolePts,HoleCenter)*n,1,3).*...
        repmat(n',length(HolePts),1);
    Nodes0(FB(:,1),:) = HolePtsProj;
end

NewNode = length(Nodes0)+1;
NewNodes = [Nodes0;HoleCenter];

NewElements = [FB,ones(length(FB),1)*NewNode];



%% Check that the normals of the newly created element are properlu oriented

U = HoleCenter-TriCenter; U = U'/norm(U);
Vctrs1 = NewNodes(NewElements(:,2),:)' - NewNodes(NewElements(:,1),:)';
Vctrs2 = NewNodes(NewElements(:,3),:)' - NewNodes(NewElements(:,1),:)';



normals = transpose(cross(Vctrs1,Vctrs2));
normals = normals./repmat(sqrt(sum(normals.^2,2)),1,3);



% Invert node ordering if the normals are inverted

if mean(normals*U)<0
    NewElements = [FB(:,1),ones(length(FB),1)*NewNode,FB(:,2)];
end



%% Write output results
NewConnectivityList = [TRin.ConnectivityList;NewElements];
TRout = triangulation(NewConnectivityList,NewNodes);

if nargout>0
    HolePtsProj;
end


end