function [ Tstring , T_str_anat , ML_Width_xp , AP_Width_xp , ProstName] = PositionProth1( SubjectCode, alpha , LongStem )
%PositionProth1 : This function place a prosthesis onto a tibia this works
%in pair with a python sript
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs :
%   SubjectCode : 'char' Code of the subject
%   alpha : 'float' valgus angle for the placement of prosthesis
%   LongStem : 'binary' Decide if a long stem should be used or not
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
tic
% addpath(strcat(pwd,'\SubFunctions'))
% addpath(strcat(pwd,'\GraphicSubFunctions'))

%% Parameters
% Prescribed Posterior slope of the implant relative to mechanical axis
beta0 = 5.0;
Zoffset_tp = 7.0; % offset of the tibial plateau plan to calculate dimension
ResectionOffset = 4.0 ;
CmtThickness = 0.0;

%% Files and folders handling
addpath(genpath(strcat(pwd,'\SubFunctions')))
addpath(genpath(strcat(pwd,'\GraphicSubFunctions')))

% Get file names and parse the mesh files to matlab
RootDir = fileparts(pwd);

%% Try to find the ".mat" file containing the mesh and associated Coordinate system

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
    [ PtMedialThirdOfTT, LegSideName, ~, ~ ] = TibialTuberosityPos(ProxTib, CS , 1);
    
    %% Save Mesh and associated CS
    save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMedialThirdOfTT','LegSideName') 
end

if LegSideName == 'R'
    Right_Knee = 1;
elseif LegSideName == 'L'
    Right_Knee = 0;
end
  
LegSide = double(2*Right_Knee - 1);


%% Perform measurements
% CS.d_tp = -CS.Origin*CS.Ztp;

% Tibial plateau size
Curve = TriPlanIntersect(ProxTib, CS.Ztp, -(CS.Origin*CS.Ztp - Zoffset_tp));
ML_Width = range(Curve.Pts*CS.Y);
AP_Width = range(Curve.Pts*CS.X);

% Varus angle
Ztp__ProjmechYZ = CS.Ztp - dot(CS.Ztp,CS.X)*CS.X;
Ztp__ProjmechYZ = Ztp__ProjmechYZ/norm(Ztp__ProjmechYZ);
Angle_Varus = rad2deg(asin(Ztp__ProjmechYZ'*CS.Y));

% Tibial Slope
Ztp__ProjmechXZ = CS.Ztp - dot(CS.Ztp,CS.Y)*CS.Y;
Ztp__ProjmechXZ = Ztp__ProjmechXZ/norm(Ztp__ProjmechXZ);
Angle_Slope = LegSide*rad2deg(asin(Ztp__ProjmechXZ'*CS.X));

% ---- Case for kinematic Alignement ---
% if a alpha angle provided is over a certain value the program
% "understand" that a kinematic alignement is prescibed 
if abs(alpha) > 40 
    alpha = -Angle_Varus;
    beta0 = abs(Angle_Slope);
    KA = 1; %Kinematic Alignment.
else
    KA = 0;
end

% Prescribe slope for the prosthesis 
beta = -LegSide*beta0;
%% Define the Cutting plan and its associated CS (U_xp, V_xp, Nxp)
% 1st define plan normal
Nxp = rot(CS.Y,beta)*rot(CS.X,alpha)*CS.Z;

% 2nd : Find minimal distance between AS and Cut plan
distMed = bsxfun(@minus,TrObjects.EpiTibASMed.Points,CS.Origin)*Nxp;
distLat = bsxfun(@minus,TrObjects.EpiTibASLat.Points,CS.Origin)*Nxp;
minDist = abs(min(min(distMed),min(distLat)));

% 3rd : Offset cut plan to get the required distance between cut plan and
% subchondral bone surface
d_xp = - ( CS.Origin - (minDist + ResectionOffset)*Nxp')*Nxp;

% Finally, get
V_xp = CS.Y - (CS.Y'*Nxp)*Nxp; V_xp= V_xp / norm(V_xp);
U_xp = cross(V_xp,Nxp);

R_xp = [U_xp V_xp Nxp];


%% Find prosthesis matching bone morphology at the resection plan
% Obtention the tibia outline at the resection plan
Curve_xp = TriPlanIntersect(ProxTib,Nxp,d_xp);
Boundary_xp = Curve_xp(1).Pts;

% Tibia dimension at prosthesis cut plan
Centroid_xp = PlanPolygonCentroid3D(Boundary_xp);
ML_Width_xp = range(Boundary_xp*V_xp);
AP_Width_xp = range(Boundary_xp*U_xp);

% Start_Point = Origin Points position of the proshtesis , [depends on prosthesis CAO,
% here the origin on CAO is located on the superior surface on the middle of the posterior edge

[ Prosthesis0, StemTip, Thickness , ProstName ] = ...
    SelectImplantSize(RootDir, ML_Width_xp, AP_Width_xp, 1, LongStem );

% Adapt to the coordinate frame and leg side [Specific of the prosthesis geometry]                 
StemTip = StemTip*[0 LegSide 0 ; 1 0 0; 0 0 -1]'; % [0 LegSide 0 ; LegSide 0 0; 0 0 -1] %[0 1 0 ; LegSide 0 0; 0 0 -1]
Prosthesis = triangulation(Prosthesis0.ConnectivityList,transpose([0 LegSide 0 ; 1 0 0; 0 0 -1]*Prosthesis0.Points'));

% Elmts2D = fixNormals( Prosthesis.Points, Prosthesis.ConnectivityList );
% Prosthesis = triangulation(Elmts2D,Prosthesis.Points);

Start_Point = Centroid_xp-U_xp'*0.36*LegSide*AP_Width_xp+0.02*V_xp'*ML_Width_xp...
    +(Thickness+CmtThickness)*Nxp'; %Prosthesis thickness +1.5 cement thickness

Oxp = Start_Point -(Thickness+CmtThickness)*Nxp';

% Move Stem Tip in CT Cooridante frame
StemTip_CT = R_xp*StemTip' + Start_Point';

CurveStemTip = TriPlanIntersect(ProxTib,Nxp,StemTip_CT);
BoundaryStemTip = CurveStemTip(1).Pts;

CDiaphysisStemTip_CT = PlanPolygonCentroid3D(BoundaryStemTip); % Center of bone at Stem Tip


%% Optimization Stem Tip Position made in the prosthesis coordinate system
Boundary_xp_inRxp = transpose(R_xp'*bsxfun(@minus,Boundary_xp,Oxp)');
Boundary_xp_inRxp = Boundary_xp_inRxp(1:7:end-1,:); 
CDiaphysisStemTip = transpose(R_xp'*(CDiaphysisStemTip_CT-Oxp)');

CurvesProsthesisTP = TriPlanIntersect(Prosthesis,[10^-6; 10^-6; 1],2); %10^-6 to avoid numerical error
Boundary_ProsthesisTP = [CurvesProsthesisTP.Pts(1:5:end-1,:) ; CurvesProsthesisTP.Pts(end,:)];


TT_on_xp = LinePlanIntersect( PtMedialThirdOfTT, CS.Z, Nxp, Oxp );
TTproj = transpose( R_xp'*(TT_on_xp' - Oxp'));

%% Optimisation of the placement of the prosthesis
    % Limit overhang
    % Orient prosthesis toward the tibial tuberosity
    % Orient prosthesis in the CS.X axis
    % if long stem, ensure that the stem tip is centered relative to the
    % diaphysis
    
     
% Geometric optimization problem
lb = [-15,-15,-15];
ub = [15,15,15];  
A = [];
b = [];
Aeq = [];
beq = [];
x0 = [0,0,0]; 
U_TT = TTproj' / norm(TTproj);

x0(3) = rad2deg(asin(U_TT(2)));

if abs(StemTip(3))>100
    f = @(x)CenteringStemTipGoal(x , Boundary_xp_inRxp , Boundary_ProsthesisTP , StemTip);  % C or CDiaphysisStemTip
    fcon = @(x)CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP , TTproj);
    x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,fcon);
else
%     f = @(x)CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP , TTproj);  % C or CDiaphysisStemTip
%     options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',500);
%     [x,fval,history] = fminunc(f,x0,options);
    [ x,fval,history ] = problem_CoverageTT( x0, Boundary_xp_inRxp, Boundary_ProsthesisTP , TTproj );
end



ProthOrig = Start_Point + x(1)*U_xp' + x(2)*V_xp';
Rp = rot(Nxp,x(3));

%% Placement Matrix To FreeCAD

PtsProsth0 = Prosthesis0.Points;
PtsProsth0(:,4) = ones(length(PtsProsth0),1);


T = zeros(4,4); T(1:3,1:3) = Rp*R_xp*[0 LegSide 0 ; LegSide 0 0; 0 0 -1]; %[0 LegSide 0 ; 1 0 0; 0 0 -1]
T(:,4)=[ProthOrig';1];


%% Plot Deformation with implanted Tibial Implant
PtsProsthEnd = transpose(T*PtsProsth0');
PtsProsthEnd(:,4)=[];

ProsthesisEnd = triangulation(Prosthesis0.ConnectivityList,PtsProsthEnd);

close all;
PlotPosOptim( ProxTib, Prosthesis0, history, Start_Point, Oxp, U_xp, V_xp, Nxp, R_xp, LegSide, d_xp, CS, PtMedialThirdOfTT, Boundary_xp, TT_on_xp, 1 )
PlotTibiaDeformation(TrObjects, ProsthesisEnd, PtMedialThirdOfTT, CS )

ProsthesisShape2 = TriPlanIntersect(Prosthesis,[10^-6; 10^-6; 1],-0.15);
[ coverage, malRotation ] = OptimOutput( x, Boundary_xp_inRxp, Boundary_ProsthesisTP, TTproj, ProsthesisShape2, LegSide );

%% Get interesting points and data to export
% Get Pts on Stem Tip planar surface
PtOnStemTip = getPtOnStemTip(ProsthesisEnd,Nxp);

% Normal at stem tip
Nst = Rp*R_xp*[0 LegSide 0 ; LegSide 0 0; 0 0 -1]*[0; sind(4.5); -cosd(4.5)];

%Nxp
%Oxp
%PtMedialThirdOfTT
% PtSemiMembranous %
% PtBiceps %
% PtSartorius % Very low forces
% PtSoleus % Insertion on fibula head and diaphysis



%% Write Results
figHandles = findobj('Type', 'figure');
figName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '.fig'];
imgName = ['Figs_' SubjectCode '_alpha' num2str(alpha) '.png'];
saveas(gcf,imgName)
savefig(figHandles,figName,'compact');
save(TmpFileName,'ProxTib','DistTib','CS','TrObjects','PtMedialThirdOfTT','LegSideName') 

outPut(end+1).Subject = SubjectCode;
outPut(end).Implant = ProstName;
outPut(end).Coverage = coverage;
outPut(end).OrientationAngle = malRotation;
outPut(end).figName = figName;
if KA
    outPut(end).alignmentType = 'KA';
else
    outPut(end).alignmentType = 'MA';
end
outPut(end).angleVV = -Angle_Varus;
outPut(end).angleTS = Angle_Slope;
outPut(end).angleAlignF = alpha;
outPut(end).angleAlignAP = beta0;

savePath = [pwd '\Output.mat'];
save(savePath,'outPut')

% 
[ CtrltyScore, Tabl ] = CentralityScore(ProxTib, Prosthesis, ProsthesisEnd, StemTip, LegSide);
writetable(Tabl,['Centrality_' SubjectCode '_alpha' num2str(alpha) '.txt'])
% 

%% Data For Dictionnary Creation
fID1=fopen(['Dict_' SubjectCode '_alpha' num2str(alpha) '.txt'],'w');
formatSpec1 = 'T, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f, %4.8f\r\n';
Tt=T';
fprintf(fID1,formatSpec1,Tt(:));
fprintf(fID1,'ZALT, %4.8f\r\n', StemTip_CT(3)+16 );
fprintf(fID1,'UDiaph, %4.8f, %4.8f, %4.8f\r\n', Zanat);
fprintf(fID1,'Xmech, %4.8f, %4.8f, %4.8f,\r\n', CS.X);
fprintf(fID1,'Ymech, %4.8f, %4.8f, %4.8f\r\n', CS.Y);
fprintf(fID1,'Zmech, %4.8f, %4.8f, %4.8f\r\n', CS.Z);
fprintf(fID1,'Nxp, %4.8f, %4.8f, %4.8f', Nxp);
fprintf(fID1,'Pt_xp, %4.8f, %4.8f, %4.8f', Oxp);
fprintf(fID1,'Pt_TT, %4.8f, %4.8f, %4.8f', PtMedialThirdOfTT);
fprintf(fID1,'Nst, %4.8f, %4.8f, %4.8f', Nst);
fprintf(fID1,'Pt_StemTip, %4.8f, %4.8f, %4.8f', PtOnStemTip);
fclose(fID1);

%% Data for FreeCAD and general data
fID3=fopen(['Output_' SubjectCode '_alpha' num2str(alpha) '.txt'],'w');
fprintf(fID3,'name= "%s" \r\n', SubjectCode );
fprintf(fID3,'ZALT= \r\n %4.8f \r\n', StemTip_CT(3)+16 );
fprintf(fID3,'Axe Diaphise : \r\n %4.8f %4.8f %4.8f  \r\n', Zanat);
fprintf(fID3,'Axe M�ca X : \r\n %4.8f %4.8f %4.8f  \r\n', CS.X);
fprintf(fID3,'Axe M�ca Y : \r\n %4.8f %4.8f %4.8f  \r\n', CS.Y);
fprintf(fID3,'Axe M�ca Z : \r\n %4.8f %4.8f %4.8f  \r\n', CS.Z);
fprintf(fID3,'Normal plateau tibial : \r\n %4.8f %4.8f %4.8f  \r\n', Nxp);
fprintf(fID3,'Angle Diaphise/Plateau tibial plan Frontal (angle varus) : \r\n %2.2f   \r\n', Angle_Varus);
fprintf(fID3,'Angle Diaphise/Plateau tibial plan Sagittal (pente tibial) : \r\n %2.2f   \r\n', Angle_Slope);

fprintf(fID3,'Centrality Scores CV : \r\n %2.2f   \r\n', CtrltyScore.CV);
fprintf(fID3,'Centrality Scores Min/Max : \r\n %2.2f   \r\n', CtrltyScore.MinMax);
fprintf(fID3,'Centrality Scores Min/Mean : \r\n %2.2f   \r\n', CtrltyScore.MinMean);

fprintf(fID3,'Partie Python pour freeCAD \r\n \r\n');
fprintf(fID3,'obj0=App.ActiveDocument.ActiveObject \r\n');
fprintf(fID3,'obj = FreeCAD.getDocument("Unnamed").getObject("Part__Feature") \r\n');
formatSpec2 = '(%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f)\r\n';
Tt=T';
fprintf(fID3,'newplace=FreeCAD.Matrix');
fprintf(fID3,formatSpec2,Tt(:));
fclose(fID3);

Tstring =sprintf(strcat('newplace=FreeCAD.Matrix',formatSpec2),T(:));

Tanat=zeros(4,4);Tanat(1:3,1:3) = CS.V*[0 LegSide 0 ; LegSide 0 0; 0 0 -1];
Tanat(:,4)=[CS.Origin';1];
Tanat = Tanat';
T_str_anat = sprintf(strcat('newplaceAnat=FreeCAD.Matrix',formatSpec2),Tanat(:));


end

