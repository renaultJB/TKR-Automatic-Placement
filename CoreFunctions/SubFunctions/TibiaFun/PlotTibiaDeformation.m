function  PlotTibiaDeformation(ProxTib, DistTib, ProsthesisEnd,  CS )
%PLOTTIBIADEFORMATION Summary of this function goes here
%   Detailed explanation goes here


%% Transform the triangulation into the tibia CS
ProxTib = triangulation(ProxTib.ConnectivityList,transpose(CS.V'*bsxfun(@minus,ProxTib.Points,CS.Origin)'));
DistTib = triangulation(DistTib.ConnectivityList,transpose(CS.V'*bsxfun(@minus,DistTib.Points,CS.Origin)'));
Prosthesis = triangulation(ProsthesisEnd.ConnectivityList,transpose(CS.V'*bsxfun(@minus,ProsthesisEnd.Points,CS.Origin)'));

%% Transform the 

Ntp = CS.V'*CS.TpCS.Z;

% Varus angle

Ztp__ProjmechYZ = Ntp - [Ntp(1); 0; 0];
Ztp__ProjmechYZ = Ztp__ProjmechYZ/norm(Ztp__ProjmechYZ);
Uvarus = [ 0; Ztp__ProjmechYZ(3); -Ztp__ProjmechYZ(2)];
Angle_Varus = rad2deg(asin(Ztp__ProjmechYZ(2)));

% Tibial Slope
Ztp__ProjmechXZ = Ntp - [0; Ntp(2); 0];
Ztp__ProjmechXZ = Ztp__ProjmechXZ/norm(Ztp__ProjmechXZ);
Uslope = [ Ztp__ProjmechXZ(3); 0; -Ztp__ProjmechXZ(1)];
% Angle_Slope = LegSide*rad2deg(asin(Ztp__ProjmechXZ(1)));



figure()
subplot(1,4,1)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
hold on
trisurf(DistTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none');
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
view(90,0)
plotCylinder( Uvarus, 1, [0 0 0], 90, 1, 'r')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(DistTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(DistTib.Points(:,3))) + 20, 1, 'k')


subplot(1,4,2)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
hold on
trisurf(DistTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none');
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( Uslope, 1, [0 0 0], 60, 1, 'r')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(DistTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(DistTib.Points(:,3))) + 20, 1, 'k')
view(0,0)

subplot(1,4,3)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.4,'edgecolor','none'); % 0.8,0.8,0.85
hold on
trisurf(DistTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.4,'edgecolor','none');
trisurf(Prosthesis,'Facecolor','g','FaceAlpha',0.4,'edgecolor','none');
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(DistTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(DistTib.Points(:,3))) + 20, 1, 'k')
view(90,0)

subplot(1,4,4)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.4,'edgecolor','none'); % 0.8,0.8,0.85
hold on
trisurf(DistTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.4,'edgecolor','none');
trisurf(Prosthesis,'Facecolor','g','FaceAlpha',0.4,'edgecolor','none');
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(DistTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(DistTib.Points(:,3))) + 20, 1, 'k')
view(0,0)

%% Figure Prox tibia


figure()
subplot(1,4,1)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
hold on
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
view(90,0)
plotCylinder( Uvarus, 1, [0 0 0], 90, 1, 'r')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(ProxTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(ProxTib.Points(:,3))) + 20, 1, 'k')


subplot(1,4,2)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
hold on
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( Uslope, 1, [0 0 0], 60, 1, 'r')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(ProxTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(ProxTib.Points(:,3))) + 20, 1, 'k')
view(0,0)

subplot(1,4,3)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.4,'edgecolor','none'); % 0.8,0.8,0.85
hold on
trisurf(Prosthesis,'Facecolor','g','FaceAlpha',0.4,'edgecolor','none');
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(ProxTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(ProxTib.Points(:,3))) + 20, 1, 'k')
view(90,0)

subplot(1,4,4)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.4,'edgecolor','none'); % 0.8,0.8,0.85
hold on
trisurf(Prosthesis,'Facecolor','g','FaceAlpha',0.4,'edgecolor','none');
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(ProxTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(ProxTib.Points(:,3))) + 20, 1, 'k')
view(0,0)


end

