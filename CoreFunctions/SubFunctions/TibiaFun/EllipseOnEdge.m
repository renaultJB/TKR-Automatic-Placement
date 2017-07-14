function [ Xel, Yel, ellipsePts, ellipse_t, EdgePtsTop ,K  ] = EllipseOnEdge( TRedge, n , d )
%EllipseOnEdge, fit an ellipse on the exterior points of the patch
%representing the medial and lateral condyles articular surfaces (AS)
%   TRedge : A triangulation object of the patch of both AS
%   A plan defined by :
%       n : normal vector ; 3x1 matrix
%       d : altitude of plan
%
% Outputs :
%   Xel : a vector (3x1 matrix) corresponding to the ellipse small axis
%   Yel : a vector (3x1 matrix) corresponding to the ellipse big axis
%   ellipsePts : Points of the ellipse (nx3 matrix)
%   ellipse_t : structure of the parameters of the fitted ellipse
%   EdgePtsTop : ????
%   K : Index of the point belonging to the convex hull.

EdgePtsID = unique(TRedge.freeBoundary);
EdgePts = TRedge.Points(EdgePtsID,:);

EdgePtsTop = EdgePts(EdgePts*n > -d - 10, : );

EdgePtsTopProj = ProjectOnPlan(EdgePtsTop,n,d);

[V,~] = eig(cov(EdgePtsTopProj));

EdgePtsTop2D = EdgePtsTopProj*V;

K = convhull(EdgePtsTop2D(:,2:3));

EdgePtsTopProjOut = EdgePtsTopProj(K,:);

ellipse_t = fit_ellipse( EdgePtsTop2D(K,2), EdgePtsTop2D(K,3) );

ellipsePts2D = [EdgePtsTop2D(1)*ones(1,length(ellipse_t.data)) ; ellipse_t.data];
Yel2D = [0 ; sin( ellipse_t.phi ) ; cos( ellipse_t.phi )];
Xel2D = [0 ; cos( ellipse_t.phi ) ; -sin( ellipse_t.phi )];



Xel = V*Xel2D;
Yel = V*Yel2D;
ellipsePts = transpose( V*ellipsePts2D );


end

