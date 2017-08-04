function [ C,Ceq ] = CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP, TT)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ThetaMax = deg2rad(10);
CutIntensity = 0.75;

x(3) = x(3)*pi/180 ;
TT(3) = [];
Boundary_xp_inRxp(:,3) = [];
Boundary_ProsthesisTP(:,3) = [];

ProsthOrig = [x(1) x(2)];

% Get the direction of the tibial tuberosity 
U_TT  = TT - ProsthOrig; U_TT = U_TT / norm(U_TT);

Theta0 = atan2(U_TT(2),U_TT(1));
if Theta0 < -pi/2   %sign(mean(Boundary_ProsthesisTP(:,1)))==-1
%     Theta = pi + atan2(U_TT(2),U_TT(1));
    Theta = pi + Theta0;
elseif Theta0 > pi/2
    Theta = Theta0 - pi;
else
    Theta = Theta0;
end

% Reduce Theta to outliers of TT orientation
Theta = sign(Theta)*(1.2198*abs(Theta)^3 - 2.3031*abs(Theta)^2 + 1.3608*abs(Theta));
Theta = 0.9*Theta;

% Identify points of the prosthesis distance to the cut plan boundary
ProsthContourT = bsxfun(@plus,Boundary_ProsthesisTP,[x(1) x(2)]) ; %Translated Contour Points

R = [cos(x(3)) -sin(x(3));sin(x(3)) cos(x(3))];
ProsthContourTR_tmp = R*bsxfun(@minus,ProsthContourT,ProsthOrig)'; % Temporary rotated and translated Contour points relative to the prosth origin R*(OPts - OProstOrig)

ProsthContourTR = bsxfun(@plus,ProsthContourTR_tmp',ProsthOrig);
d = p_poly_dist(ProsthContourTR(:,1), ProsthContourTR(:,2), Boundary_xp_inRxp(:,1), Boundary_xp_inRxp(:,2));

% Criterium 1 : No overhang of the implant
C1 = 30 * mean(exp(d+0.5).^2);  % Cost function , d+0.5 for penalty if the prosthesis is too close of the edges

% Criterium 2 : Implant oriented towards the Tibial tuberosity
deltaTheta = rad2deg(x(3))-rad2deg(Theta);
C2 = (1 - exp(-(deltaTheta)^2/100))*(deltaTheta)^2;

% Criterium 3 : orientation vs cut plan;
C3 = 200*abs(x(3));

% Overall, sum of the three criteria
C = C1 + C2 + C3;

Ceq = -1;




% pause(0.01)
% U_prosth = sign(mean(Boundary_ProsthesisTP(:,1)))*[cos(x(3)) x(3)];
% 
% figure(8)
% plot(ProsthContourTR(:,1),ProsthContourTR(:,2))
% hold on
% plot(Boundary_xp_inRxp(:,1),Boundary_xp_inRxp(:,2),'.-k')
% plot(TT(1),TT(2),'rs')
% plot([TT(1);x(1)],[TT(2);x(2)],'r.-')
% plot([x(1);x(1)+20*U_prosth(1)],[x(2);x(2)+20*U_prosth(2)],'b*-')
% axis equal
% hold off


end

