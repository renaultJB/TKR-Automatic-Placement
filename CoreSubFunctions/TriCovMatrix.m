function [ CovM ] = TriCovMatrix( TR )
%Compute the covariance matrix of a triangulation by taking the triangle
% centers and weighting them by their corresponding area

Ppt = TriMesh2DProperties(TR);

W = Ppt.Area;
WTot = Ppt.TotalArea;

CovM = zeros(3);

Pts = TR.incenter;

for i=1:3
    for j=1:3
        CovM(i,j) = sum((Pts(:,i)-mean(Pts(:,i))).*(Pts(:,j)-mean(Pts(:,j))).*W)/WTot;
    end
end



end

