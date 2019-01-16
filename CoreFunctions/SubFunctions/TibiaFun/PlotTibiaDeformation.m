function  PlotTibiaDeformation(ProxTib, DistTib, ProsthesisEnd,  CS )
%PLOTTIBIADEFORMATION Summary of this function goes here
%   Detailed explanation goes here


%% Transform the triangulation into the tibia CS
ProxTib = triangulation(ProxTib.ConnectivityList,transpose(CS.V'*bsxfun(@minus,ProxTib.Points,CS.Origin)'));
DistTib = triangulation(DistTib.ConnectivityList,transpose(CS.V'*bsxfun(@minus,DistTib.Points,CS.Origin)'));
Prosthesis = triangulation(ProsthesisEnd.ConnectivityList,transpose(CS.V'*bsxfun(@minus,ProsthesisEnd.Points,CS.Origin)'));
EpiTibArtMed = triangulation(CS.Morph.EpiTibArtMed.ConnectivityList,...
    transpose(CS.V'*bsxfun(@minus,CS.Morph.EpiTibArtMed.Points,CS.Origin)'));
%% Transform the 

Ntp = CS.V'*CS.TpCS.Z; %

[NtpMed,d] = LS_Plan( EpiTibArtMed.Points );



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
% Medial Tibial Slope
ZtpMed__ProjmechXZ = NtpMed - [0; NtpMed(2); 0];
ZtpMed__ProjmechXZ = ZtpMed__ProjmechXZ/norm(ZtpMed__ProjmechXZ);
UslopeMed = [ ZtpMed__ProjmechXZ(3); 0; -ZtpMed__ProjmechXZ(1)];

AltEnd = -[0,0,1]*CS.V'*(CS.Origin*CS.Z-CS.AltEndDiaph)*CS.Z;
AltStart = min(ProxTib.Points*[0;0;1]);


Alt = AltStart+4 : AltEnd;
FrontEdgePts=ones(length(Alt),3);


i=0;
for alt=Alt
    i=i+1;
    Curve = TriPlanIntersect(ProxTib, [0;0;1], alt);
    [~,Imax] = min(Curve(1).Pts(:,1));
%     Imax
%     Curve(1).Pts(Imax,:)
    PostEdgePts(i,1:3) = Curve(1).Pts(Imax,:);
end


PostEdgePts = ProjectOnPlan(PostEdgePts,[0;1;0],zeros(1,3));

[V,~] = eig(cov(PostEdgePts));

UEdge = V(:,3); UEdge=sign(UEdge(3))*UEdge;

CenterEdgePost = mean(PostEdgePts);





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


% subplot(1,4,2)  
% trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
% hold on
% trisurf(DistTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none');
% axis equal
% grid off
% axis off
% light('Position',[50 50 200],'Style','local')
% light('Position',[-50 -50 -200],'Style','local')
% light('Position',[-50 100 -200],'Style','local')
% plotCylinder( Uslope, 1, [0 0 0], 80, 1, 'r')
% plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(DistTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(DistTib.Points(:,3))) + 20, 1, 'k')
% plotCylinder( UEdge, 0.75, CenterEdgePost + 25*UEdge', 180, 1, 'b')
% view(180,0)

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
plotCylinder( UslopeMed, 1, mean(EpiTibArtMed.Points), 80, 1, 'r')
plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(DistTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(DistTib.Points(:,3))) + 20, 1, 'k')
plotCylinder( UEdge, 0.75, CenterEdgePost + 25*UEdge', 180, 1, 'b')
view(180,0)

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


% subplot(1,4,2)
% trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
% hold on
% axis equal
% grid off
% axis off
% light('Position',[50 50 200],'Style','local')
% light('Position',[-50 -50 -200],'Style','local')
% light('Position',[-50 100 -200],'Style','local')
% plotCylinder( Uslope, 1, [0 0 0], 80, 1, 'r')
% % plotCylinder( [0; 0; -1], 0.75, [0, 0, 0.5*min(ProxTib.Points(:,3)) + 0.5*max(ProxTib.Points(:,3))], abs(min(ProxTib.Points(:,3))) + 20, 1, 'k')
% plotCylinder( UEdge, 0.75, CenterEdgePost + 25*UEdge' , 180, 1, 'b')
% view(180,0)

subplot(1,4,2)
trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
hold on
axis equal
grid off
axis off
light('Position',[50 50 200],'Style','local')
light('Position',[-50 -50 -200],'Style','local')
light('Position',[-50 100 -200],'Style','local')
plotCylinder( UslopeMed, 1, mean(EpiTibArtMed.Points), 80, 1, 'r')
plotCylinder( UEdge, 0.75, CenterEdgePost + 25*UEdge', 180, 1, 'b')
view(180,0)

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

