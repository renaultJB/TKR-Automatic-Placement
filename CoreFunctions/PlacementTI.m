function [ Tstring , T_str_anat , ML_Width_xp , AP_Width_xp , ProstName, alpha] = PlacementTI( SubjectCode, alpha , implantType , LongStem, beta )
%PlacementTI : This function place a tibial implant onto a tibia this works
%in pair with a python sript
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs :
%   SubjectCode : 'char' Code of the subject
%   alpha : 'float' valgus angle for the placement of prosthesis
%   LongStem : 'binary' Decide if a long stem should be used or not
%   implantType : 'char' Type of implant for the database we have
%
% Outputs :
%   Tstring
%   T_str_anat
%   ML_Width
%   Ap_Width
%   ProstName
%
% Lexic :
%   TP : Tibial Plateau
%   0 (at the end of variables) : initial value
%   CS : List of coordinate system of the tibia
%   Rct : CT Scan coordinate system
%   Rxp : Coordinate system of the tibial cut plan
%   xp : cutting plane
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Parameters
% Prescribed Posterior slope of the implant relative to mechanical axis
if nargin == 4
    beta = 7;
end
Zoffset_tp = 10; % offset of the tibial plateau plan to calculate dimension
ResectionOffset = 3. ; % 2mm + 1mm saw blade thickness
CmtThickness = 0.;
PhysioTTAangle = 12.5;

%% Files and folders handling
addpath(genpath(strcat(pwd,'\SubFunctions')))
addpath(genpath(strcat(pwd,'\GraphicSubFunctions')))

% Get file names and parse the mesh files to matlab
RootDir = fileparts(pwd);

%% Try to find the ".mat" file containing the mesh and already computed
%  anatomical Coordinate system

TmpFileName = strcat(RootDir,'\tempFiles\',SubjectCode,'.mat');
TmpFileExist = exist(TmpFileName,'file');

if TmpFileExist ~= 0
    load(TmpFileName)
else
    ProxTibMeshFile = strcat(RootDir,'\Tibia_',SubjectCode,'.msh');
    DistTibMeshFile = strcat(RootDir,'\DistTibia_',SubjectCode,'.msh');
    
    %% Read mesh files of the proximal an distal tibia
    [ProxTib,DistTib] = ReadMesh( ProxTibMeshFile, DistTibMeshFile );
    
    %% Construct the coordinates system of the tibia
    [ CSs, TrObjects ] = RTibiaFun( ProxTib , DistTib);
    CS = CSs.PIAASL;

    %% Find the tibial tuberosity of the Tibia, this also permits the identification of the legside  
    [ PtMedialThirdOfTT, LegSideName, ~, ~, PtMiddleOfTT ] = TibialTuberosityPos(ProxTib, CS , 1);
       
    %% Save Mesh and associated CS
    save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMedialThirdOfTT','PtMiddleOfTT','LegSideName') 
end

%% Get leg side code 
if LegSideName == 'R'
    Right_Knee = 1;
elseif LegSideName == 'L'
    Right_Knee = 0;
end
LegSide = double(2*Right_Knee - 1);

%% Projection matrices on tibia anatomic coordinate system
% Frontal Plan projection
Af = [CS.Y,CS.Z];
CS.Pfront = Af*inv(Af'*Af)*Af';

% Saggital plan projection
As = [CS.Z,CS.X];
CS.Psag = As*inv(As'*As)*As';

% Axial (transverse) plan projection
At = [CS.X,CS.Y];
CS.Paxial = At*inv(At'*At)*At';

%% Perform measurements
% "Tibial Varus" angle (MMPTA)
Ztp__ProjmechYZ = normalizeV(CS.Pfront*CS.Ztp);
Angle_Varus = rad2deg(asin(Ztp__ProjmechYZ'*CS.Y));
mMPTA = 90 - Angle_Varus;

% Tibial Slope
Angle_Slope = -LegSide*rad2deg(asin(CS.Ztp'*CS.X));




CS.Anat.Alpha = -Angle_Varus;
CS.Anat.Beta = Angle_Slope;

CS.nTP =   - LegSide * sin(deg2rad(CS.Anat.Beta))*CS.X ...
    - cos(deg2rad(CS.Anat.Beta))*sin(deg2rad(CS.Anat.Alpha))*CS.Y ...
    + cos(deg2rad(CS.Anat.Beta))*cos(deg2rad(CS.Anat.Alpha))*CS.Z;
% ---- Case for kinematic Alignement ---
% if a alpha angle provided is over a certain value the program
% "understand" that a kinematic alignement is prescibed
if abs(alpha) > 30
    alpha = -Angle_Varus;
    beta = Angle_Slope;
    KA = 1; %Kinematic Alignment.
    Alignment = 'Kine';
else
    KA = 0;
    Alignment = 'Mech';
end

%% Define the Cutting plan and its associated CS (U_xp, V_xp, Nxp)
% 1st define plan normal, 
%   For right knee beta indicate the anterior rotation of the
%   cut plan, since the clinical convention indicate positive
%   tibial slope for posterior slope, we introduce the LegSide factor 
Nxp =   - LegSide * sin(deg2rad(beta))*CS.X ...
        - cos(deg2rad(beta))*sin(deg2rad(alpha))*CS.Y ...
        + cos(deg2rad(beta))*cos(deg2rad(alpha))*CS.Z;


    
% 2nd : Find minimal distance between AS and Cut plan
distMed = bsxfun(@minus,TrObjects.EpiTibASMed.Points,CS.Origin)*Nxp;
distLat = bsxfun(@minus,TrObjects.EpiTibASLat.Points,CS.Origin)*Nxp;
minDist = min(min(distMed),min(distLat));


% 3rd : Offset cut plan to get the required distance between cut plan and
% subchondral bone surface, Oc a point on the cut plan
Oxp = CS.Origin + (minDist - ResectionOffset)*Nxp';
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

%% Plots
[ PtMedialThirdOfTT, LegSideName, ~, ~, PtMiddleOfTT ] = TibialTuberosityPos(ProxTib, CS , 1);
hold on; plotArrow(Nxp,1.5,CS.Origin,30,1,'r')
plotArrow(CS.Z,1,CS.Origin,30,1,'k')
pl3t(Boundary_xp,'g-')

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
exitflag = -1000 ;
reduceSize = 0 ;
tries = 0 ;
while exitflag < -1 && tries < 3
    switch implantType
        case {'Nexgen','nexgen',1}
            prosth_type = 1;
            [ Prosthesis0, StemTip, Thickness , ProstName ] = ...
                SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, prosth_type, LongStem, LegSideName, reduceSize  );
            
            %         TI_speTransfo = [0 LegSide 0 ; 1 0 0; 0 0 -1]; % [0 LegSide 0 ; LegSide 0 0; 0 0 -1]
            TI_speTransfo = [0 LegSide 0 ; LegSide 0 0; 0 0 -1] ;
            % Normal at stem Tip plan
            Nst_0 = [0; sind(4.5); -cosd(4.5)];
        case {'Persona','persona',3}
            prosth_type = 3;
            [ Prosthesis0, StemTip, Thickness , ProstName ] = ...
                SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, prosth_type, LongStem, LegSideName, reduceSize );
            
            TI_speTransfo = [0 -LegSide 0 ; LegSide 0 0; 0 0 1];
            
            
        case {'GC6','GrandChallenge','GC',4}
            prosth_type = 4;
            [ Prosthesis0, StemTip, Thickness , ProstName ] = ...
                SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, prosth_type, LongStem, LegSideName, reduceSize );
            
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
    
    
    %
    % figure(50)
    % trisurf(Prosthesis)
    % hold on
    % pl3t(CurvesProsthesisTP.Pts,'r*-')
    
    % Simplify Implant
    Boundary_ProsthesisTP = [CurvesProsthesisTP.Pts(1:6:end-1,:) ; CurvesProsthesisTP.Pts(end,:)];
    
    
    % figure(51)
    % pl3t(CurvesProsthesisTP.Pts,'r-')
    % axis equal
    % hold on
    % pl3t(Boundary_xp_inRxp,'g-')
    % plotDot([0,0,0],'r',0.5)
    % plotDot(PtMiddleOfTT_inRxp,'g',0.5)
    
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
    
    % figure(10)
    % pl3t(Boundary_xp_inRxp,'k-')
    % hold on
    % axis equal
    % pl3t(Boundary_ProsthesisTP,'r-')
    
    [ x,fval,exitflag,output,history ] = optimC_PlacementTI_xp( x0, A,b,Aeq,beq,lb,ub, Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA );
    % [ x,fval,history ] = optimUC_PlacementTI_xp( x0, Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA );
    
    if exitflag < -1 
        reduceSize = 1 ;
    end
    tries = tries + 1 ;
    
end


ProthOrig = Start_Point + x(1)*U_xp' + x(2)*V_xp';
Rp = rot(Nxp,x(3));

% Get Delta Theta
U_it = [cosd(x(3)) ; sind(x(3)) ; 0];
U_t = normalizeV(CS.Paxial*(R_xp*U_it));
theta_it = rad2deg(acos(CS.Y'*U_t));
deltaTheta = theta_TTA - theta_it; % Rotational error in degree 18° is the physiological value

%
gamma = deg2rad(x(3));
R = [cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0; 0 0 1];
ProsthContourTR_tmp = R*Boundary_ProsthesisTP';
        % 2nd translate origin
ProsthContourTR_optim = bsxfun(@plus,ProsthContourTR_tmp',[x(1) x(2) 0]);

pl3t(ProsthContourTR_optim,'b-')

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

%% Placement Matrix of the implant
PtsProsth0 = Prosthesis0.Points;
PtsProsth0(:,4) = ones(length(PtsProsth0),1);

T = zeros(4,4); T(1:3,1:3) = Rp*R_xp*TI_speTransfo; %[0 LegSide 0 ; 1 0 0; 0 0 -1]
T(:,4)=[ProthOrig';1];

PtsProsthEnd = transpose(T*PtsProsth0');
PtsProsthEnd(:,4)=[];

ProsthesisEnd = triangulation(Prosthesis0.ConnectivityList,PtsProsthEnd);

close all;
%% Plot Deformation with implanted Tibial Implant
% 
% PlotPosOptim( ProxTib, Prosthesis0, history, Start_Point, Oxp, U_xp, V_xp, R_xp, CS, PtMiddleOfTT, Boundary_xp, TI_speTransfo)
% PlotTibiaDeformation(TrObjects, ProsthesisEnd, PtMiddleOfTT, CS )


%% Coverage and malrotation relativ to AP Axis
areaTibia = polyarea(Boundary_xp_inRxp(:,1),Boundary_xp_inRxp(:,2));
ProsthContourTR = TriPlanIntersect(ProsthesisEnd,Nxp,Oxp+0.5*Nxp');
ProsthContourTR = ProsthContourTR(1).Pts;
areaImplant = polyarea(ProsthContourTR(:,1),ProsthContourTR(:,2));
coverage = 100*areaImplant/areaTibia;

malRotation = x(3);
%% Get interesting points and data to export
% Get Pts on Stem Tip planar surface
PtOnStemTip = getPtOnStemTip(ProsthesisEnd,Nxp);

% Normal at stem tip
Nst = T(1:3,1:3)*Nst_0;

% Get central point of medial and lateral condyles
PtLat_0 = CS.CenterKnee - 0.275*CS.ML_Width_xp*CS.Y';
IdPt = ProxTib.nearestNeighbor(PtLat_0);
PtLat = ProxTib.Points(IdPt,:);

PtMed_0 = CS.CenterKnee + 0.275*CS.ML_Width_xp*CS.Y';
IdPt = ProxTib.nearestNeighbor(PtMed_0);
PtMed = ProxTib.Points(IdPt,:);

% Posterior point for semimembranosous
Pt_SM_0 = CS.CenterKnee ...
          + 0.325*CS.ML_Width_xp*CS.Y' ...
          - 0.5*CS.AP_Width_xp*LegSide*CS.X' ...
          - 0.125*(CS.ML_Width_xp+CS.AP_Width_xp)*CS.Z' ;
IdPt = ProxTib.nearestNeighbor(Pt_SM_0);
Pt_SM = ProxTib.Points(IdPt,:);

% Pt Fibula
Pt_Fib_0 = CS.CenterKnee ...
          - 0.5*CS.ML_Width_xp*CS.Y' ...
          - 0.4*CS.AP_Width_xp*LegSide*CS.X' ...
          - 0.22*(CS.ML_Width_xp+CS.AP_Width_xp)*CS.Z' ;
IdPt = ProxTib.nearestNeighbor(Pt_Fib_0);
Pt_Fib = ProxTib.Points(IdPt,:);

% Pt SemiTendinosous
Pt_ST_0 = PtMiddleOfTT ...
          + 0.2*CS.ML_Width_xp*CS.Y' ...
          - 0.1*CS.AP_Width_xp*LegSide*CS.X'  ;
IdPt = ProxTib.nearestNeighbor(Pt_ST_0);
Pt_ST = ProxTib.Points(IdPt,:);

% Pt Soleus
Pt_Sol_0 = PtMiddleOfTT ...
          - 0.05*CS.ML_Width_xp*CS.Y' ...
          - 0.9*CS.AP_Width_xp*LegSide*CS.X'...
          - 0.8*CS.AP_Width_xp*CS.Z'  ;
IdPt = ProxTib.nearestNeighbor(Pt_Sol_0);
Pt_Sol = ProxTib.Points(IdPt,:);

% Pt Tibialis Anterior
Pt_TA_0 = PtMiddleOfTT ...
          - 0.15*CS.ML_Width_xp*CS.Y' ...
          - 0.5*CS.AP_Width_xp*LegSide*CS.X'...
          - 0.6*CS.AP_Width_xp*CS.Z';
IdPt = ProxTib.nearestNeighbor(Pt_TA_0);
Pt_TA = ProxTib.Points(IdPt,:);

% Pt Tibialis Posterior
Pt_TP_0 = PtMiddleOfTT ...
          - 0.1*CS.ML_Width_xp*CS.Y' ...
          - 0.9*CS.AP_Width_xp*LegSide*CS.X'...
          - 1.0*CS.AP_Width_xp*CS.Z'   ;
IdPt = ProxTib.nearestNeighbor(Pt_TP_0);
Pt_TP = ProxTib.Points(IdPt,:);

%% Write ouptut files
[ CtrltyScore, Tabl ] = CentralityScore(ProxTib, Prosthesis, ProsthesisEnd, StemTip);
writetable(Tabl,['Centrality_' SubjectCode '_alpha' num2str(alpha) '.txt'])


%% Data For Dictionnary Creation
fID1=fopen(['Dict_' SubjectCode '_alpha' num2str(alpha) '.txt'],'w');
fprintf(fID1,'ImpltType, "%s"\r\n', ProstName );
fprintf(fID1,'LegSide, "%s"\r\n', LegSideName );
fprintf(fID1,'Alignment, "%s"\r\n', Alignment );
formatSpec1 = 'T, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f\r\n';
Tt=T';
fprintf(fID1,formatSpec1,Tt(:));
fprintf(fID1,'ZALT, %4.1f\r\n', StemTip_CT(3)+16 );
fprintf(fID1,'alpha, %4.2f \r\n', alpha );
fprintf(fID1,'beta, %4.2f \r\n', beta );% Posterior slope of implantation relative to mechanical axis
fprintf(fID1,'mMPTA, %4.2f \r\n', mMPTA );
fprintf(fID1,'TPslope, %4.2f \r\n', Angle_Slope ); % Posterior slope of tibial plateau relative to mechanical axis
fprintf(fID1,'UDiaph, %4.8f, %4.8f, %4.8f\r\n', CS.Z0);
fprintf(fID1,'O_knee, %4.8f, %4.8f, %4.8f\r\n', CS.CenterKnee);
fprintf(fID1,'Xmech, %4.8f, %4.8f, %4.8f\r\n', CS.X);
fprintf(fID1,'Ymech, %4.8f, %4.8f, %4.8f\r\n', CS.Y);
fprintf(fID1,'Zmech, %4.8f, %4.8f, %4.8f\r\n', CS.Z);
fprintf(fID1,'Ztp, %4.8f, %4.8f, %4.8f\r\n', CS.Ztp);
fprintf(fID1,'Nxp, %4.8f, %4.8f, %4.8f\r\n', Nxp);
fprintf(fID1,'Pt_xp, %4.8f, %4.8f, %4.8f\r\n', Oxp);
fprintf(fID1,'Pt_LCC, %4.8f, %4.8f, %4.8f\r\n', PtLat);
fprintf(fID1,'Pt_MCC, %4.8f, %4.8f, %4.8f\r\n', PtMed);
fprintf(fID1,'Pt_MTTT, %4.8f, %4.8f, %4.8f\r\n', PtMedialThirdOfTT);
fprintf(fID1,'Pt_TT, %4.8f, %4.8f, %4.8f\r\n', PtMiddleOfTT); 
fprintf(fID1,'Pt_Fib, %4.8f, %4.8f, %4.8f\r\n', Pt_Fib); % Tibia fibular articular surface center
fprintf(fID1,'Pt_SM, %4.8f, %4.8f, %4.8f\r\n', Pt_SM); % Semimembranosous tibial insertion
fprintf(fID1,'Pt_ST, %4.8f, %4.8f, %4.8f\r\n', Pt_ST); % Semitendinosous tibial insertion
fprintf(fID1,'Pt_Sol, %4.8f, %4.8f, %4.8f\r\n', Pt_Sol); % Soleus tibial insertion
fprintf(fID1,'Pt_TA, %4.8f, %4.8f, %4.8f\r\n', Pt_TA); % Tibialis Anterior tibial insertion
fprintf(fID1,'Pt_TP, %4.8f, %4.8f, %4.8f\r\n', Pt_TP); % Tibialis Posterior tibial insertion
fprintf(fID1,'Nst, %4.8f, %4.8f, %4.8f\r\n', Nst);
fprintf(fID1,'Pt_StemTip, %4.3f, %4.3f, %4.3f\r\n', PtOnStemTip);
fclose(fID1);

%% Data for FreeCAD and general data
fID3=fopen(['Output_' SubjectCode '_alpha' num2str(alpha) '.txt'],'w');
fprintf(fID3,'name= "%s" \r\n', SubjectCode );
fprintf(fID3,'ZALT= \r\n %4.8f \r\n', StemTip_CT(3)+16 );
fprintf(fID3,'alpha= \r\n %4.8f \r\n', alpha );
fprintf(fID3,'Axe Diaphise : \r\n %4.8f %4.8f %4.8f  \r\n', CS.Z0);
fprintf(fID3,'Axe Meca X : \r\n %4.8f %4.8f %4.8f  \r\n', CS.X);
fprintf(fID3,'Axe Meca Y : \r\n %4.8f %4.8f %4.8f  \r\n', CS.Y);
fprintf(fID3,'Axe Meca Z : \r\n %4.8f %4.8f %4.8f  \r\n', CS.Z);
fprintf(fID3,'Normal plateau tibial : \r\n %4.8f %4.8f %4.8f  \r\n', Nxp);
fprintf(fID3,'Angle Diaphise/Plateau tibial plan Frontal (angle varus) : \r\n %2.2f   \r\n', Angle_Varus);
fprintf(fID3,'Angle Diaphise/Plateau tibial plan Sagittal (pente tibiale) : \r\n %2.2f   \r\n', Angle_Slope);
fprintf(fID3,'Coverage (percent) : \r\n %2.2f   \r\n', coverage);
fprintf(fID3,'Rotation (degrees) : \r\n %2.2f  \r\n', malRotation);
fprintf(fID3,'Centrality Scores CV : \r\n %2.2f   \r\n', CtrltyScore.CV);
fprintf(fID3,'Centrality Scores Min_Max : \r\n %2.2f   \r\n', CtrltyScore.MinMax);
fprintf(fID3,'Centrality Scores Min_Mean : \r\n %2.2f   \r\n', CtrltyScore.MinMean);

fprintf(fID3,'Partie Python pour freeCAD \r\n \r\n');
fprintf(fID3,'obj0=App.ActiveDocument.ActiveObject \r\n');
fprintf(fID3,'obj = FreeCAD.getDocument("Unnamed").getObject("Part__Feature") \r\n');
formatSpec2 = '(%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f)\r\n';
Tt=T';
fprintf(fID3,'newplace=FreeCAD.Matrix');
fprintf(fID3,formatSpec2,Tt(:));
fclose(fID3);

Tstring =sprintf(strcat('newplace=FreeCAD.Matrix',formatSpec2),Tt(:));

Tanat=zeros(4,4);Tanat(1:3,1:3) = CS.V*[0 LegSide 0 ; LegSide 0 0; 0 0 -1];
Tanat(:,4)=[CS.Origin';1];
Tanat = Tanat';
T_str_anat = sprintf(strcat('newplaceAnat=FreeCAD.Matrix',formatSpec2),Tanat(:));


%% Write Results
% figHandles = findobj('Type', 'figure');
% figName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '.fig'];
% imgName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '.png'];
% saveas(gcf,imgName)
% savefig(figHandles,figName,'compact');
% save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMiddleOfTT','LegSideName') 

% outPut(end+1).Subject = SubjectCode;
% outPut(end).Implant = ProstName;
% outPut(end).Coverage = coverage;
% outPut(end).OrientationAngle = malRotation;
% outPut(end).figName = figName;
% if KA
%     outPut(end).alignmentType = 'KA';
% else
%     outPut(end).alignmentType = 'MA';
% end
% outPut(end).angleVV = -Angle_Varus;
% outPut(end).angleTS = Angle_Slope;
% outPut(end).angleAlignF = alpha;
% outPut(end).angleAlignAP = beta;
% 
% savePath = [pwd '\Output.mat'];
% save(savePath,'outPut')


% figure()
% trisurf(Prosthesis)
% hold on
% pl3t(CurvesProsthesisTP.Pts,'r*-')
% Boundary_ProsthesisTP = [CurvesProsthesisTP.Pts(1:5:end-1,:) ; CurvesProsthesisTP.Pts(end,:)];
% 
% 

%% Export Pos File for control meshing in GMSH
Mat3D_final = writeGMSHPosFile(ProxTib,ProsthesisEnd,2.2,1.7,2.7,Nxp,Oxp);
%save('Mat3D_final.mat','Mat3D_final','ProxTib','ProsthesisEnd','Nxp','Oxp')

%% Close all

close all

end

