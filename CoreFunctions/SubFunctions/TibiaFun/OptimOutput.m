function [ coverage, malRotation ] = OptimOutput( x, Boundary_xp_inRxp, Boundary_ProsthesisTP, TTproj, ProsthesisShape2, LegSide )
%OPTIMOUTPUT Compute Coverage and malrotation angle

Boundary_xp_inRxp = Boundary_xp_inRxp(:,1:2);
Boundary_ProsthesisTP = Boundary_ProsthesisTP(:,1:2);
TTproj(3) = [];
x(3) = deg2rad(x(3));

%% Malrotation
ProsthOrig = [x(1) x(2)];
ProsthContourT = bsxfun(@plus,Boundary_ProsthesisTP,ProsthOrig);
R = [cos(x(3)) -sin(x(3));sin(x(3)) cos(x(3))];
ProsthContourTR_tmp = R*bsxfun(@minus,ProsthContourT,ProsthOrig)';
ProsthContourTR = bsxfun(@plus,ProsthContourTR_tmp',ProsthOrig);


for i=1:length(ProsthesisShape2)
    ProsthesisShape2(i).Pts = ProsthesisShape2(i).Pts(:,1:2);
    
    ProsthesisShape2(i).Pts = bsxfun(@plus,ProsthesisShape2(i).Pts,ProsthOrig);
    R = [cos(x(3)) -sin(x(3));sin(x(3)) cos(x(3))];
    ProsthesisShape2(i).Pts = R*bsxfun(@minus,ProsthesisShape2(i).Pts,ProsthOrig)';
    ProsthesisShape2(i).Pts = bsxfun(@plus,ProsthesisShape2(i).Pts',ProsthOrig);
    
end

LineTT = [x(1) x(2); TTproj ];
ImplantAxis = LegSide*[0 ,0; cos(x(3)) , sin(x(3))];
ImplantAxis = 55*ImplantAxis + [ProsthOrig; ProsthOrig];

U_TT = diff(LineTT)/norm(diff(LineTT));
U_Implant = diff(ImplantAxis)/norm(diff(ImplantAxis));


figure()
plot(Boundary_xp_inRxp(:,1),Boundary_xp_inRxp(:,2),'g-','lineWidth',2)
hold on 
axis equal
plot(ProsthContourTR(:,1),ProsthContourTR(:,2),'r-','lineWidth',1.5)
for i=1:length(ProsthesisShape2)
    plot(ProsthesisShape2(i).Pts(:,1),ProsthesisShape2(i).Pts(:,2),'r-','lineWidth',1)
end
plot(TTproj(1),TTproj(2),'k*')

plot(LineTT(:,1),LineTT(:,2),'g-','lineWidth',1)
plot(LineTT(:,1),LineTT(:,2),'k.-','lineWidth',1)

plot(ImplantAxis(:,1),ImplantAxis(:,2),'r-','lineWidth',1)
plot(ImplantAxis(:,1),ImplantAxis(:,2),'k.-','lineWidth',1)


malRotation = acosd(U_TT*U_Implant');

%% Coverage
areaTibia = polyarea(Boundary_xp_inRxp(:,1),Boundary_xp_inRxp(:,2));
areaImplant = polyarea(ProsthContourTR(:,1),ProsthContourTR(:,2));

coverage = areaImplant/areaTibia;
end

