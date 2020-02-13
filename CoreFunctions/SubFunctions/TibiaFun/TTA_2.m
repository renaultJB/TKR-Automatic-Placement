
% Alternative way to identify the TTA using projection along the mechanical
% axis, unfinished, paused project...

AltitudesFract = 0.8: 0.1 : 1.2;

D_Epi = CS.Origin*CS.Z - CS.AltStartEpiph;
D_Dia = CS.Origin*CS.Z - CS.AltEndDiaph;


ProxTib_inRT = triangulation(ProxTib.ConnectivityList,...
    transpose(CS.V'*bsxfun(@minus,ProxTib.Points,CS.Origin)'));

[ PtMedialThirdOfTT, LegSideName, ~, ~, PtMiddleOfTT ] = TibialTuberosityPos(ProxTib_inRT, CS , 1);

U_mech = [10^-6; 10^-6; 1];
Curve = TriPlanIntersect(ProxTib_inRT,U_mech,1.1*D_Epi);
CSPts0 = Curve.Pts;

Curve = TriPlanIntersect(ProxTib_inRT,U_mech,D_Dia);
CSPts0 = Curve.Pts;

Curve = TriPlanIntersect(ProxTib_inRT,U_mech,0.25*D_Epi);
CSPts1 = Curve.Pts;


figure()
trisurf(ProxTib_inRT,'Facecolor',[0.65    0.65    0.6290],...
    'FaceAlpha',1,'edgecolor','none'); % 0.8,0.8,0.85
hold on
axis equal
light('Position', 300*CS.Y + 200*CS.X,'Style','local')
light('Position', 200*CS.Y - 200*CS.X,'Style','local')
light('Position', 50*CS.Y + 50*CS.X - 500*CS.Z,'Style','local')

hold on
grid off
lighting gouraud
pl3t(CSPts0,'g-')
pl3t(CSPts1,'r-')


figure()
pl3t(CSPts0,'g-')
hold on;
pl3t(CSPts1,'r-')


% plotDot( PtsMedThird, 'r', 1.25 )
% plotDot( PtMedialThirdOfTT, 'g', 2.5 )
% plotDot( PtsTT, 'b', 1.25 )