function [ C,Ceq ] = CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP, TT)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

U_prosth = sign(mean(Boundary_ProsthesisTP(:,1)))*[cos(x(3)) x(3)];

ProsthContourT = bsxfun(@plus,Boundary_ProsthesisTP,[x(1) x(2)]) ; %Translated Contour Points

R = [cos(x(3)) -sin(x(3));sin(x(3)) cos(x(3))];
ProsthContourTR_tmp = R*bsxfun(@minus,ProsthContourT,ProsthOrig)'; % Temporary rotated and translated Contour points relative to the prosth origin R*(OPts - OProstOrig)

ProsthContourTR = bsxfun(@plus,ProsthContourTR_tmp',ProsthOrig);

d = p_poly_dist(ProsthContourTR(:,1), ProsthContourTR(:,2), Boundary_xp_inRxp(:,1), Boundary_xp_inRxp(:,2));

bProsth = poly2mask(ProsthContourTR(:,1)-min(Boundary_xp_inRxp(:,1))+50, ProsthContourTR(:,2)-min(Boundary_xp_inRxp(:,2))+50,150,150);
bCut = poly2mask(Boundary_xp_inRxp(:,1)-min(Boundary_xp_inRxp(:,1))+50, Boundary_xp_inRxp(:,2)-min(Boundary_xp_inRxp(:,2))+50,150,150);

bOutsideArea = bProsth - bCut.*bProsth;

% imshow(bOutsideArea);

AreaOut = sum(sum(bOutsideArea));

% Criterium 1 : No overhang of the implant
D = d(d>0);

if isempty(D)
    D=0;
end

C1 = 25 * mean(exp(d).^2);  %+ exp((x(3)+0.1)^4/5)-1; %exp(d).^2-1+0.5*d % Cost function , d+2 for penalty if the prosthesis is too close of the edges

% Criterium 2 : Implant oriented towards the Tibial tuberosity
% C2 = 25 * (sqrt(1-(U_TT*U_prosth')^2))

C2 = (rad2deg(x(3))-rad2deg(Theta))^2;

rad2deg(x(3))
rad2deg(Theta)
% Criterium 3 : Area outside the tibial cut section;
% C3 = AreaOut

% Overall, sum of the two criteria
C = C1 + C2; % + C3;


% Add a part on the orientation relative to the TTA : 
% (x(3)-xTTA)^6/100    exp((x+5.7)^2/5)-1


Ceq = -1;

rad2deg(x(3))

pause(0.1)

figure(1)
plot(ProsthContourTR(:,1),ProsthContourTR(:,2))
hold on
plot(Boundary_xp_inRxp(:,1),Boundary_xp_inRxp(:,2),'.-k')
plot(TT(1),TT(2),'rs')
plot([TT(1);x(1)],[TT(2);x(2)],'r.-')
plot([x(1);x(1)+20*U_prosth(1)],[x(2);x(2)+20*U_prosth(2)],'b*-')
axis equal
hold off


end

