% retrieve Data
T_op = [0.1521022024836078 0.9871624325776751 0.0487365541145966 -89.3514978457129700 ;
    -0.9844438324435510 0.1557024120114873 -0.0814070000532018 20.3279791537453178 ;
    -0.0879503312301400 -0.0355962161069904 0.9954886481700262 -30.4056042081937967 ;
    0.0000000000000000 0.0000000000000000 0.0000000000000000 1.0000000000000000 ];

O_xp_op = [-89.9264 23.3014 -38.0441];

O_xp_op = [-96.3609 17.4345 -38.7190];
Nxp_op = T_op(1:3,3);
O_xp_op = CS.Origin - ((CS.Origin - O_xp_op)*Nxp_op)*Nxp_op';
Thickness = 8.39;

V_xp_op = normalizeV( CS.Y - (CS.Y'*Nxp_op)*Nxp_op );
U_xp_op = cross(V_xp_op,Nxp_op);
R_xp_op = [U_xp_op V_xp_op Nxp_op];



beta_op = asind(-LegSide*Nxp_op'*CS.X);
alpha_op = asind(-Nxp_op'*CS.Y/cos(deg2rad(beta_op)));

% Get gamma
R_IT_op = T_op(1:3,1:3);
Rp = R_xp_op'*R_IT_op*TI_speTransfo';
gamma = asind(Rp(2,1));
x(3) = gamma;

% Get Delta Theta
U_it = [cos(x(3)) ; sin(x(3)) ; 0];
U_t = normalizeV(CS.Paxial*(R_xp*U_it));
theta_it = rad2deg(acos(CS.Y'*U_t));
deltaTheta = theta_TTA - theta_it; % Rotational error in degree 18° is the physiological value


% Get Start Point
ProthOrig = T_op(1:3,4)';

Start_Point_op = O_xp_op + (Thickness+CmtThickness)*Nxp_op';
x(1) = (ProthOrig - Start_Point_op)*U_xp_op;
x(2) = (ProthOrig - Start_Point_op)*V_xp_op;


%% Plot Results
history = x;

Curve_xp = TriPlanIntersect(ProxTib,Nxp_op,O_xp_op);
Boundary_xp = Curve_xp(1).Pts;

PlotPosOptim2D( x, Prosthesis, Boundary_xp, ProxTib, GS, GS_TTA,...
    PtMiddleOfTT, O_xp_op, PtMedialThirdOfTT, R_xp_op, CS,  history(end,:)  )


figHandles = findobj('Type', 'figure');
figName = ['Figs_' SubjectCode '_alpha' num2str(alpha_op) '_beta' num2str(beta_op) '.fig'];
imgName = ['Figs_' SubjectCode '_alpha' num2str(alpha_op) '_beta' num2str(beta_op) '.png'];
saveas(gcf,imgName)
savefig(figHandles,figName,'compact');
save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMiddleOfTT','LegSideName')



PlotPosOptim( ProxTib, Prosthesis0, history, Start_Point_op, O_xp_op, U_xp_op, V_xp_op, R_xp_op, CS, PtMiddleOfTT, Boundary_xp, TI_speTransfo)


close all

%% Get best placement in resection plane
% Copy normal code but change alpha beta and Oxp
beta = beta_op
alpha = alpha_op

%% Normal code
Nxp =   - LegSide * sin(deg2rad(beta))*CS.X ...
        - cos(deg2rad(beta))*sin(deg2rad(alpha))*CS.Y ...
        + cos(deg2rad(beta))*cos(deg2rad(alpha))*CS.Z;


    
% 2nd : Find minimal distance between AS and Cut plan
distMed = bsxfun(@minus,TrObjects.EpiTibASMed.Points,CS.Origin)*Nxp;
distLat = bsxfun(@minus,TrObjects.EpiTibASLat.Points,CS.Origin)*Nxp;
minDist = min(min(distMed),min(distLat));


% 3rd : Offset cut plan to get the required distance between cut plan and
% subchondral bone surface, Oc a point on the cut plan
% Oxp = CS.Origin + (minDist - ResectionOffset)*Nxp';
Oxp = O_xp_op
d_xp = - Oxp*Nxp ;

% Finally, get
V_xp = normalizeV( CS.Y - (CS.Y'*Nxp)*Nxp );
U_xp = cross(V_xp,Nxp);
R_xp = [U_xp V_xp Nxp];

%% Get the geometry of the tibia at the resection plan
% Obtention the tibia outline at the resection plan
Curve_xp = TriPlanIntersect(ProxTib,Nxp,d_xp);
Boundary_xp = Curve_xp(1).Pts;

% Tibia dimension at prosthesis cut plan
ML_Width_xp = range(Boundary_xp*V_xp);
AP_Width_xp = range(Boundary_xp*U_xp);

CS.ML_Width_xp = ML_Width_xp;
CS.AP_Width_xp = AP_Width_xp;

%% Get the TTA orientation with Berger et al. method
% Get the tibia contour at the proximal measurement
% 5% under the sqrt of the area at the cut level
O_pm = Oxp - 0.05*sqrt(pi/4*ML_Width_xp*AP_Width_xp)*CS.Z';
Curve_pm = TriPlanIntersect(ProxTib,CS.Z,O_pm);

Boundary_pm_inRt = transpose(CS.V'*Curve_pm(1).Pts');

% Get boundary above the cutting plan
ellipse_t = fit_ellipse( Boundary_pm_inRt(:,1), ...
                Boundary_pm_inRt(:,2));
GS_inRt = mean(ellipse_t.data,2);

GS = CS.V*[GS_inRt(1);GS_inRt(2);mean(Boundary_pm_inRt(:,3))];

Curve_ttam = TriPlanIntersect(ProxTib,CS.Z,PtMiddleOfTT);

GS_TTA = GS - ((GS'-PtMiddleOfTT)*CS.Z)*CS.Z;
U_TTA = normalizeV(PtMiddleOfTT' - GS_TTA);
theta_TTA = acosd(CS.Y'*U_TTA) - PhysioTTAangle;
theta_TTA = mod(180 + LegSide*theta_TTA , 180);



GS_MTTTA = GS - ((GS'-PtMedialThirdOfTT)*CS.Z)*CS.Z;
U_MTTTA = normalizeV(PtMedialThirdOfTT' - GS_TTA);
theta_MTTTA = acosd(CS.Y'*U_MTTTA);
theta_MTTTA = mod(180 + LegSide*theta_MTTTA , 180);

% theta_TTA = theta_MTTTA;

%% Find prosthesis matching bone morphology at the resection plan
% Start_Point = Origin Points position of the proshtesis , [depends on prosthesis CAO,
% here the origin on CAO is located on the superior surface on the middle of the posterior edge

%% Switch between specifities of implants
switch implantType
    case {'Nexgen','nexgen',1}
        prosth_type = 1;
        [ Prosthesis0, StemTip, Thickness , ProstName ] = ...
            SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, prosth_type, LongStem );
        
        %         TI_speTransfo = [0 LegSide 0 ; 1 0 0; 0 0 -1]; % [0 LegSide 0 ; LegSide 0 0; 0 0 -1]
        TI_speTransfo = [0 LegSide 0 ; LegSide 0 0; 0 0 -1] ;
    case {'Persona','persona',3}
        prosth_type = 3;
        [ Prosthesis0, StemTip, Thickness , ProstName ] = ...
            SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, prosth_type, LongStem, LegSideName );
        
        TI_speTransfo = [0 -LegSide 0 ; LegSide 0 0; 0 0 1];
        
        
    case {'GC6','GrandChallenge','GC',4}
        prosth_type = 4;
        [ Prosthesis0, StemTip, Thickness , ProstName ] = ...
            SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, prosth_type, LongStem, LegSideName );
        
        TI_speTransfo = [LegSide 0 0 ; 0 LegSide 0 ; 0 0 1];
        
        
    otherwise
        warning('Type of implants not implemented wet')
end


StemTip = StemTip*TI_speTransfo';
Prosthesis = triangulation(Prosthesis0.ConnectivityList,...
    transpose(TI_speTransfo*Prosthesis0.Points'));

Start_Point = Oxp + (Thickness+CmtThickness)*Nxp';

% Move Stem Tip in CT Coordinate frame
StemTip_CT = R_xp*StemTip' + Start_Point';

CurveStemTip = TriPlanIntersect(ProxTib,Nxp,StemTip_CT);
BoundaryStemTip = CurveStemTip(1).Pts;

% Center of bone at Stem Tip :
CDiaphysisStemTip_CT = PlanPolygonCentroid3D(BoundaryStemTip); 

%% Optimization of the implant position in the resection plan coordinate system
Boundary_xp_inRxp = transpose(R_xp'*bsxfun(@minus,Boundary_xp,Oxp)');
Boundary_xp_inRxp = Boundary_xp_inRxp(1:8:end-1,:); % Reduce Boundary Matriw weight
CDiaphysisStemTip = transpose(R_xp'*(CDiaphysisStemTip_CT-Oxp)');
PtMiddleOfTT_inRxp = transpose(R_xp'*(PtMiddleOfTT-Oxp)');

CurvesProsthesisTP = TriPlanIntersect(Prosthesis,[10^-6; 10^-6; 1],2.5); %10^-6 to avoid numerical error

if length(CurvesProsthesisTP) > 1
    AreaMax=0;
    for i = 1 : length(CurvesProsthesisTP)
        A = polyarea(CurvesProsthesisTP(i).Pts(:,1),CurvesProsthesisTP(i).Pts(:,2));
        if A > AreaMax
            AreaMax = A;
            iMax = i;
        end
    end
    CurvesProsthesisTP = CurvesProsthesisTP(iMax);
end

CurvesProsthesisTP.Pts = bsxfun(@plus,CurvesProsthesisTP.Pts, [0,0,2.5]);



figure(50)
trisurf(Prosthesis)
hold on
pl3t(CurvesProsthesisTP.Pts,'r*-')

% Simplify Implant 
Boundary_ProsthesisTP = [CurvesProsthesisTP.Pts(1:6:end-1,:) ; CurvesProsthesisTP.Pts(end,:)];


figure(51)
pl3t(CurvesProsthesisTP.Pts,'r-')
axis equal
hold on
pl3t(Boundary_xp_inRxp,'g-')
plotDot([0,0,0],'r',0.5)
plotDot(PtMiddleOfTT_inRxp,'g',0.5)

%% Optimization of the placement of the prosthesis
    % Limit overhang
    % Orient prosthesis toward the tibial tuberosity
    % Orient prosthesis in the CS.X axis
    % if long stem, ensure that the stem tip is centered relative to the
    % diaphysis

     
% Geometric optimization problem
lb = [-30,-20,-20];
ub = [30,20,20];  
A = [1,1,0;-1,-1,0;];
b = [100,100];
Aeq = [];
beq = [];
x0 = [0,0,0]; 


% initial value guess
O_it0 = PlanPolygonCentroid3D(Boundary_ProsthesisTP); 
x0(1) = - O_it0(1);
x0(2) = - O_it0(2);
x0(3) = theta_TTA - 18 - 90;
x0(3) = 0;

[ x,fval,history ] = optimC_PlacementTI_xp( x0, A,b,Aeq,beq,lb,ub, Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA );   
% [ x,fval,history ] = optimUC_PlacementTI_xp( x0, Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA );

ProthOrig = Start_Point + x(1)*U_xp' + x(2)*V_xp';
Rp = rot(Nxp,x(3));

% Get Delta Theta
U_it = [cos(deg2rad(x(3))) ; sin(deg2rad(x(3))) ; 0];
U_t = normalizeV(CS.Paxial*(R_xp*U_it));
theta_it = rad2deg(acos(CS.Y'*U_t));
deltaTheta = theta_TTA - theta_it; % Rotational error in degree 18° is the physiological value



gamma = deg2rad(x(3));
R = [cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0; 0 0 1];
ProsthContourTR_tmp = R*Boundary_ProsthesisTP';
        % 2nd translate origin
ProsthContourTR_optim = bsxfun(@plus,ProsthContourTR_tmp',[x(1) x(2) 0]);


%% Plot 2D placement with TTA in Rt
% Initial State
% close all
PlotPosOptim2D( x, Prosthesis, Boundary_xp, ProxTib, GS, GS_TTA,...
    PtMiddleOfTT, Oxp, PtMedialThirdOfTT, R_xp, CS, history(1,:) )
figHandles = findobj('Type', 'figure');
% figName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '.fig'];
imgName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '_beta' num2str(beta) '_init.png'];
saveas(gcf,imgName)
% savefig(figHandles,figName,'compact');
save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMiddleOfTT','LegSideName')


% Post Optimiasation
close all
PlotPosOptim2D( x, Prosthesis, Boundary_xp, ProxTib, GS, GS_TTA,...
    PtMiddleOfTT, Oxp, PtMedialThirdOfTT, R_xp, CS,  history(end,:)  )


figHandles = findobj('Type', 'figure');
figName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '_beta' num2str(beta) '.fig'];
imgName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '_beta' num2str(beta) '.png'];
saveas(gcf,imgName)
savefig(figHandles,figName,'compact');
save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMiddleOfTT','LegSideName')

PlotPosOptim( ProxTib, Prosthesis0, history, Start_Point, Oxp, U_xp, V_xp, R_xp, CS, PtMiddleOfTT, Boundary_xp, TI_speTransfo)


