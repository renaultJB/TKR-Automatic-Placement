function [ Pt_intersection ] = LinePlanIntersect( oLine, uLine, nPlan, oPlan )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Inputs :
%   - oLine : a [1x3] matrix of a point on the line
%   - uLine : a [3x1] Matrix of the directing vector of the line
%   - nPlan : a [3x1] Matrix of the vector normal to the plan
%   - oPlan : a [1x3] matrix of a point on the plan
%
% Output :
%   - Pt_intersection : a [1x3] matrix, coordinates of intersection point
%
%
%
% -------------------------------------------------------------------------

if  uLine'*nPlan == 0
    error('No or non-finite intersection, the line direction and the plan normal are perpendiculars')
end

t = ( (oPlan-oLine)*nPlan ) / ( uLine'*nPlan );

Pt_intersection = oLine + t * uLine';


% Matrix oriented alternative

% A = -eye(4);
% A(16) = 0;
% A(4,1:3)= nPlan';
% A(1:3,4) = uLine;
% 
% B = - oLine';
% B(4) = oPlan*nPlan; 
% 
% X = A\B;
% 
% Pt_intersection = X(1:3);

end