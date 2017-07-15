function [ Tstring , T_str_anat , PC_ML_Width , PC_AP_Width , ProstName] = PositionProth1( SubjectCode, alpha, Right_Knee , LongStem , Generate_Pos )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
close all
% Lexic
% TP : Tibial Plateau
% 0 (at the end) : initial value
% Rct : CT Scan coordinate system
% Rtb : Tibial eigen vector coordinate system
% 2D : Refer to the surface mesh of the object
% xp : cutting plane
tic
addpath(strcat(pwd,'\SubFunctions'))

%% Parameters
LegSide = double(2*Right_Knee - 1); % 1 for right knee, -1 for left

r=5;
Taille_Elmt_Base = 5; % en mm
Taille_Elmt_Raf = 0.3; % (a peu pr�s) en mm


%% Get file names and parse the mesh files to matlab

cd ../
ProxTibMeshFile = strcat(pwd,'\Tibia_',SubjectCode,'.msh');
DistTibMeshFile = strcat(pwd,'\DistTibia_',SubjectCode,'.msh');
cd ./CoreFunctions

[ ProxTib, DistTib ] = ReadCheckMesh( ProxTibMeshFile, DistTibMeshFile );

[ CS ] = TibiaCS( ProxTib , DistTib);

[ PtMedialThirdOfTT, LegSide, PtsMedThird, PtsTT ] = TibialTuberosityPos(ProxTib, CS , 1);

%% Get a pseudo latero-medial vector obtained from the the medial anterior plan on the tibial diaphysis
% Create a 1st approximation of the medio-lateral vector & Post Ant Vector
% : XLatMed_0, YPostAnt_0
Pts_TT = Slice_pop(Pts2D,VZ,Centroid_Diaphysis,3);
sqrt(sum(bsxfun(@minus,Pts_TT,Centroid_Diaphysis).^2,2));
[~,I] = max(sqrt(sum(bsxfun(@minus,Pts_TT,Centroid_Diaphysis).^2,2)));
V_TT = Pts_TT(I,:)-Centroid_Diaphysis; V_TT = V_TT'/norm(V_TT);

if dot(V_TT,V_TP(:,1))>0
    XLatMed_0 = - V_TP(:,1);
else
    XLatMed_0 = V_TP(:,1);
end

if dot(V_TT,V_TP(:,2))<0
    YPostAnt_0 = - V_TP(:,2);
else
    YPostAnt_0 = V_TP(:,2);
end

%% Identify the TT
% Calculate the metaphysis
 Pts2D_meta = Pts2D(bsxfun(@minus,Pts2D,Centroid_End_Diaphysis')*VZ>-10 &...
     bsxfun(@minus,Pts2D,Centroid_TP')*VZ<0 , : );
 
% Tibial Tuberosity LandMark    
TT_lm = TT_position( Pts2D_meta,V_TT,VZ );




%Filter  edge points too far from condyle plan
Centroid_TP_onPlan = (Centroid_TP -(Centroid_TP' - [0 0 -dc/nc(3)])*nc*nc)';
Pedge_TP = Pedge(bsxfun(@minus,Pedge,Centroid_TP_onPlan)*nc<3 & ...
    bsxfun(@minus,Pedge,Centroid_TP_onPlan)*nc > -7,:);

% Assign Pedge to the lateral or medial condyles
PedgeMed = Pedge_TP(bsxfun(@minus,Pedge_TP,Centroid_TP_onPlan)*XLatMed_2>TP_Central/2,:);
PedgeLat = Pedge_TP(bsxfun(@minus,Pedge_TP,Centroid_TP_onPlan)*XLatMed_2<-TP_Central/2,:);

%% Calculate new TP Slice for inertial axis
%Find extremities of Tibial plateau edge in pseudo-frontal plan
[~,I0Med]  = max(bsxfun(@minus,PedgeMed,Centroid_TP')*XLatMed_2);
[~,I0Med(2)]  = min(bsxfun(@minus,PedgeMed,Centroid_TP')*XLatMed_2+...
    10*abs(bsxfun(@minus,PedgeMed,Centroid_TP')*YPostAnt_2).^2);
Pt_TP_limit_Med = PedgeMed(I0Med,:);
IMed = rangesearch(PedgeMed,Pt_TP_limit_Med,2.5);
Pts_TP_limit_Med = PedgeMed([IMed{1} IMed{2}],:);

[~,I0Lat]  = min(bsxfun(@minus,PedgeLat,Centroid_TP')*XLatMed_2);
[~,I0Lat(2)]  = max(bsxfun(@minus,PedgeLat,Centroid_TP')*XLatMed_2-...
    10*abs(bsxfun(@minus,PedgeLat,Centroid_TP')*YPostAnt_2).^2);
Pt_TP_limit_Lat = PedgeLat(I0Lat,:);
ILat = rangesearch(PedgeLat,Pt_TP_limit_Lat,2.5);
Pts_TP_limit_Lat = PedgeLat([ILat{1} ILat{2}],:);

ML_Width = norm(Pt_TP_limit_Med(1,:)-Pt_TP_limit_Lat(1,:));

Seeds = [Pts_TP_limit_Med;Pts_TP_limit_Lat];

% Calculate Pts on condyles

PtsCondyleMed = EdgeCondyle( PcondyleMed , PedgeMed, Pts2D , Seeds );
PtsCondyleLat = EdgeCondyle( PcondyleLat , PedgeLat, Pts2D , Seeds );

% Merge Condyles obtained from the 2 methods
PcondyleMed = unique([PtsCondyleMed;PcondyleMed],'rows');
PcondyleLat = unique([PtsCondyleLat;PcondyleLat],'rows');

Pcondyle = [PcondyleLat;PcondyleMed];
[nc,dc] = LS_Plan(Pcondyle);

% Centroids
CMed = mean( PcondyleMed );
CLat = mean( PcondyleLat );
Centroid_TP_onPlan = (Centroid_TP -(Centroid_TP' - [0 0 -dc/nc(3)])*nc*nc)';

% Get the axe
XLatMed_3  = CMed - CLat;
XLatMed_3  = XLatMed_3' / norm(XLatMed_3);


% PcondyleLat = PtsCondyleLat;
%% Compute the anatomical based on tibial plateau morphology
Zanat = nc;
Yanat = cross(Zanat,XLatMed_3)/norm(cross(Zanat,XLatMed_3));
Xanat = cross(Yanat,Zanat);
if YPostAnt_0'*Yanat<0
    Yanat = - Yanat;
    warning('Yanat was flipped and is Ranat is no longer a direct orthonormed coordinate system')
end
Vanat=[Xanat Yanat Zanat];

%% Compute Mechanical tibial Coordinate System


% Calculate subject varus preOP
Zanat_ProjmechXZ = Zanat - dot(Zanat,Ymech)*Ymech; Zanat_ProjmechXZ=Zanat_ProjmechXZ/norm(Zanat_ProjmechXZ);
Angle_Varus = rad2deg(atan2( sign(dot(cross(Zmech,Zanat_ProjmechXZ),Ymech))*norm(cross(Zmech,Zanat_ProjmechXZ)),dot(Zmech,Zanat_ProjmechXZ))) ;


% Determine offset of Condyle plan to cutting plan
distML = abs(atan(deg2rad(alpha - Angle_Varus)))*ML_Width/2 ;


%% Obtention of Normal of the cutting plane (Nrml_xp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = -LegSide*5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Zprost = rot(Xmech,beta)*rot(Ymech,alpha)*Zmech;

if YPostAnt_0'*Ymech<0
    Ymech = - Ymech;
    warning('Ymech was flipped and is Ranat is no longer a direct orthonormed coordinate system')
end

% Tibial plateau slope
Zanat_ProjmechYZ = Zanat - dot(Zanat,Xmech)*Xmech; Zanat_ProjmechYZ=Zanat_ProjmechYZ/norm(Zanat_ProjmechYZ);
Angle_Slope = rad2deg(atan2( norm(cross(Zmech,Zanat_ProjmechYZ)),dot(Zmech,Zanat_ProjmechYZ))); % Supposed always positive

distAP = abs(atan(deg2rad(Angle_Slope-abs(beta))))*ML_Width/4;

Yprost = cross(Zprost , XLatMed_3) / norm(cross(Zprost , XLatMed_3));
Xprost = cross(Yprost,Zprost);

if YPostAnt_0'*Yprost<0
    Yprost = - Yprost;
    warning('Yprost was flipped and is Rptost is no longer a direct orthonormed coordinate system')
end

Vprost = [Xprost Yprost Zprost]; % Coordinate base of prosthesis relative to CT coordinate System

% Adjustment of depth according to alpha angle

d_xp = - ( Centroid_TP_onPlan - (max(distML,distAP)+5)*Zprost')*Zprost;

Points_PlanCut = Slice_pop(Points,Zprost,d_xp,0.5);

% Tibia dimension at prosthesis cut plan
Centroid_pc = mean(Points_PlanCut);
PC_ML_Width = range((bsxfun(@minus,Points_PlanCut,Centroid_pc)*Xprost));
PC_AP_Width = range((bsxfun(@minus,Points_PlanCut,Centroid_pc)*Yprost));

% Start_Point = Origin Points position of the proshtesis , [depends on prosthesis CAO,
% here the origin on CAO is located on the superior surface on the middle of the posterior edge

CWD = cd; %Current Working Directory
[ Pts_prosthesis_init, StemTip, Thickness , ProstName ] = ...
    SelectImplantSize(CWD, PC_ML_Width, PC_AP_Width, 1, LongStem );


Start_Point = Centroid_pc-Yprost'*0.36*PC_AP_Width+0.02*LegSide*Xprost'*PC_ML_Width...
    +(Thickness+1.5)*Zprost'; %Prosthesis thickness +1.5 cement thickness      +LegSide*alpha*Vpc(:,1)'



% Move Stem Tip in CT Cooridante frame
StemTip_CT = Vprost*[LegSide 0 0; 0 1 0; 0 0 -1]*StemTip' + Start_Point';
StemTip = transpose([-1 0 0; 0 1 0; 0 0 -1]*StemTip');
PtsStemTip = Slice_pop(Points,Zmech,StemTip_CT,1.5); % Points of diaphysis stem tips
CDiaphysisStemTip_CT = mean(PtsStemTip) ; % Center of bone at Stem Tip


%% Optimization Stem Tip Position made in the prosthesis coordinate system

% Shape of the tibial plateau contour at the cut plan
TibialCutPts = [Slice_pop(Pts2D_BU,Zprost,d_xp,1) ; Slice_pop(Points,Zprost,d_xp,1)];
% TibialCutPts = Slice_pop(Points,Zprost,d_xp,1);
TibialCutPts = transpose( Vprost'*TibialCutPts' ) ; 
CDiaphysisStemTip = transpose( Vprost'*CDiaphysisStemTip_CT' ) ;

kTC = boundary(TibialCutPts(:,1),TibialCutPts(:,2),1);
% kTC = convhull(TibialCutPts(:,1),TibialCutPts(:,2));
kTC = kTC(1:5:end); %Downsampling boundary for performance

TibialCutBoundaries = [TibialCutPts(kTC,1) TibialCutPts(kTC,2)];

TibialCutBoundaries = [TibialCutPts(kTC,1) TibialCutPts(kTC,2)];

Pts_prosthesis0 = transpose([-1 0 0; 0 1 0; 0 0 -1]*Pts_prosthesis_init');
ProsthOrig0 = transpose( Vprost'*Start_Point');
Pts_prosthesis = bsxfun(@plus , Pts_prosthesis0 , ProsthOrig0);

kP = convhull(Pts_prosthesis(:,1),Pts_prosthesis(:,2));
ProsthContour = [Pts_prosthesis(kP,1) Pts_prosthesis(kP,2)];

ProsthOrig0_proj = Start_Point * Vprost;


TTproj = TT_lm' * Vprost; TTproj(3) = ProsthOrig0_proj(3);

% Geometric optimization problem
lb = [-15,-15,-15];
ub = [15,15,15];
A = [];
b = [];
Aeq = [];
beq = [];
x0 = [0,0,0]; 

U_TT  = TTproj - ProsthOrig0; U_TT = U_TT / norm(U_TT);

% x0(3) = -sign(U_TT(1))*rad2deg(acos(U_TT(2))) + (-5 + 10*rand(1));
x0(3) = -sign(U_TT(1))*rad2deg(acos(U_TT(2)))/2;


if StemTip(3)>100
    f = @(x)CenteringStemTipGoal(x , ProsthOrig0 , CDiaphysisStemTip , StemTip);  % C or CDiaphysisStemTip
    fcon = @(x)CoverageCost(x, TibialCutBoundaries, ProsthContour , ProsthOrig0_proj);
    x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,fcon);
else
    f = @(x)CoverageCost2(x, TibialCutBoundaries, ProsthContour , ProsthOrig0_proj, TTproj);  % C or CDiaphysisStemTip
    options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',500);
    [x,fval] = fminunc(f,x0,options);
end
fval

if fval > 150;
    [ Pts_prosthesis_init, StemTip, Thickness , ProstName ] = ...
    SelectImplantSize(CWD, PC_ML_Width, PC_AP_Width, 1, LongStem, 1 );


    Start_Point = Centroid_pc-Yprost'*0.36*PC_AP_Width+0.02*LegSide*Xprost'*PC_ML_Width...
        +(Thickness+1.5)*Zprost'; %Prosthesis thickness +1.5 cement thickness      +LegSide*alpha*Vpc(:,1)'



    % Move Stem Tip in CT Cooridante frame
    StemTip_CT = Vprost*[LegSide 0 0; 0 1 0; 0 0 -1]*StemTip' + Start_Point';
    StemTip = transpose([-1 0 0; 0 1 0; 0 0 -1]*StemTip');
    PtsStemTip = Slice_pop(Points,Zmech,StemTip_CT,1.5); % Points of diaphysis stem tips
    CDiaphysisStemTip_CT = mean(PtsStemTip) ; % Center of bone at Stem Tip


    %% Optimization Stem Tip Position made in the prosthesis coordinate system
    % Shape of the tibial plateau contour at the cut plan
    TibialCutPts = [Slice_pop(Pts2D_BU,Zprost,d_xp,1) ; Slice_pop(Points,Zprost,d_xp,1)];
    % TibialCutPts = Slice_pop(Points,Zprost,d_xp,1);
    TibialCutPts = transpose( Vprost'*TibialCutPts' ) ; 
    CDiaphysisStemTip = transpose( Vprost'*CDiaphysisStemTip_CT' ) ;


    kTC = boundary(TibialCutPts(:,1),TibialCutPts(:,2),1);
    % kTC = convhull(TibialCutPts(:,1),TibialCutPts(:,2));
    kTC = kTC(1:5:end); %Downsampling boundary for performance

TibialCutBoundaries = [TibialCutPts(kTC,1) TibialCutPts(kTC,2)];

    TibialCutBoundaries = [TibialCutPts(kTC,1) TibialCutPts(kTC,2)];

    Pts_prosthesis0 = transpose([-1 0 0; 0 1 0; 0 0 -1]*Pts_prosthesis_init');
    ProsthOrig0 = transpose( Vprost'*Start_Point');
    Pts_prosthesis = bsxfun(@plus , Pts_prosthesis0 , ProsthOrig0);

    kP = convhull(Pts_prosthesis(:,1),Pts_prosthesis(:,2));
    ProsthContour = [Pts_prosthesis(kP,1) Pts_prosthesis(kP,2)];

    ProsthOrig0_proj = Start_Point * Vprost;


    TTproj = TT_lm' * Vprost; TTproj(3) = ProsthOrig0_proj(3);

    % Geometric optimization problem
    lb = [-15,-15,-15];
    ub = [15,15,15];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    x0 = [0,0,0]; 

    if StemTip(3)>100
        f = @(x)CenteringStemTipGoal(x , ProsthOrig0 , CDiaphysisStemTip , StemTip);  % C or CDiaphysisStemTip
        fcon = @(x)CoverageCost(x, TibialCutBoundaries, ProsthContour , ProsthOrig0_proj);
        x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,fcon);
    else
        f = @(x)CoverageCost2(x, TibialCutBoundaries, ProsthContour , ProsthOrig0_proj, TTproj);  % C or CDiaphysisStemTip
        options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',500);
        [x,fval] = fminunc(f,x0,options);
    end
    
end


% x(1) = -LegSide*x(1);
x(3) = -LegSide*x(3);



ProthOrig = Start_Point + x(1)*Xprost' + x(2)*Yprost';
R = rot(Zprost,x(3));

%% Placement Matrix To FreeCAD

T=zeros(4,4); T(1:3,1:3)=R*Vprost*[LegSide 0 0; 0 1 0; 0 0 -1];
T(:,4)=[ProthOrig';1];

fID3=fopen(['Output_' name '_alpha' num2str(alpha) '.txt'],'w');
fprintf(fID3,'name= "%s" \r\n', name);
fprintf(fID3,'Axe Diaphise : \r\n %4.8f %4.8f %4.8f  \r\n', Zanat);
fprintf(fID3,'Axe M�ca X : \r\n %4.8f %4.8f %4.8f  \r\n', Xmech);
fprintf(fID3,'Axe M�ca Y : \r\n %4.8f %4.8f %4.8f  \r\n', Ymech);
fprintf(fID3,'Axe M�ca Z : \r\n %4.8f %4.8f %4.8f  \r\n', Zmech);
fprintf(fID3,'Normal plateau tibial : \r\n %4.8f %4.8f %4.8f  \r\n', Zprost);
fprintf(fID3,'Angle Diaphise/Plateau tibial plan Frontal (angle varus) : \r\n %2.2f   \r\n', Angle_Varus);
fprintf(fID3,'Angle Diaphise/Plateau tibial plan Sagittal (pente tibial) : \r\n %2.2f   \r\n', Angle_Slope);
fprintf(fID3,'Partie Python pour freeCAD \r\n \r\n');
fprintf(fID3,'obj0=App.ActiveDocument.ActiveObject \r\n');
fprintf(fID3,'obj = FreeCAD.getDocument("Unnamed").getObject("Part__Feature") \r\n');
% formatSpec2 = '((%4.8f,%4.8f,%4.8f,%4.8f),(%4.8f,%4.8f,%4.8f,%4.8f),(%4.8f,%4.8f,%4.8f,%4.8f),(%4.8f,%4.8f,%4.8f,%4.8f))\r\n';
formatSpec2 = '(%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f,%4.8f)\r\n';
Tt=T';
fprintf(fID3,'newplace=FreeCAD.Matrix');
fprintf(fID3,formatSpec2,Tt(:));
fclose(fID3);

Tstring =sprintf(strcat('newplace=FreeCAD.Matrix',formatSpec2),Tt(:));

Tanat=zeros(4,4);Tanat(1:3,1:3) = Vanat*[LegSide 0 0; 0 1 0; 0 0 -1];
Tanat(:,4)=[ProthOrig';1];
Tanat = Tanat';
T_str_anat = sprintf(strcat('newplaceAnat=FreeCAD.Matrix',formatSpec2),Tanat(:));

%% Writting of a .pos file (post processing field de gmsh)
% Use it as a background mesh, to control mesh size of the peri-prosthetic
% zone
if Generate_Pos==1

    %% surface points to 3D image representing the surface
    Pts_prosthesis = Pts_prosthesis_init;

    Tmp_Pts = ones(size(Pts_prosthesis,1),size(Pts_prosthesis,2)+1);
    Tmp_Pts(1:size(Pts_prosthesis,1),1:size(Pts_prosthesis,2)) = Pts_prosthesis;
    Tmp_Pts=(T*Tmp_Pts')';
    Pts_prosthesis=Tmp_Pts(1:size(Pts_prosthesis,1),1:size(Pts_prosthesis,2));

    %% surface points to 3D image representing the surface
    Pts2D = Pts2D_BU; % Recall original file read because we suppress some in the previous processes

    Xb_init=Pts2D(:,1); Yb_init=Pts2D(:,2); Zb_init=Pts2D(:,3);
    Xp_init=Pts_prosthesis(:,1); Yp_init=Pts_prosthesis(:,2); Zp_init=Pts_prosthesis(:,3);

    % Vertex of the Bounding box composed of the tibia and the prosthesis
    Xmin=min([min(Xp_init) min(Xb_init)]);
    Ymin=min([min(Yp_init) min(Yb_init)]);
    Zmin=min([min(Zp_init) min(Zb_init)]);

    Xmax=max([max(Xp_init) max(Xb_init)]);
    Ymax=max([max(Yp_init) max(Yb_init)]);
    Zmax=max([max(Zp_init) max(Zb_init)]);

    Xb_trslt=Xb_init-Xmin;
    Yb_trslt=Yb_init-Ymin;
    Zb_trslt=Zb_init-Zmin;

    Xp_trslt=Xp_init-Xmin;
    Yp_trslt=Yp_init-Ymin;
    Zp_trslt=Zp_init-Zmin;

    Mat3D=ones(ceil(Xmax-Xmin)+20,ceil(Ymax-Ymin)+20,ceil(Zmax-Zmin)+20).*(1/Taille_Elmt_Base);
    Length_Z = ceil(Zmax-Zmin)+20;


    % Weighting of surface
    clearvars k z Z
    for k=1:length(Pts2D)
        z = round(Zb_trslt(k))+10;
        Mat3D(round(Xb_trslt(k))+10,round(Yb_trslt(k))+10,round(Zb_trslt(k))+10)=...
            1/sqrt(Taille_Elmt_Raf)*(atan(-(0.65*Length_Z-z)/(Length_Z/16))/(pi/2)+1.01);
        Z(k) = z;
    end

    for l=1:length(Pts_prosthesis)
        Mat3D(round(Xp_trslt(l))+10,round(Yp_trslt(l))+10,round(Zp_trslt(l))+10)=2.5/Taille_Elmt_Raf;
    end


    Mat3D_blurred = imgaussfilt3(Mat3D,4.05);

    Mat3D_final=Mat3D_blurred.^-1;


    
    
    fID2=fopen(['outmsh_a' num2str(alpha) '.pos'],'w');
    fprintf(fID2,'View "background mesh" {\r\n');
    
    %format du fichier .pos : SH pour scalar hexahedron
    formatSpec = 'SH(%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f){%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f};\r\n';
    for x=1:2:ceil(Xmax-Xmin)+18
        x/(ceil(range(Xp_trslt))+18)*100 % id�e de l'avancement de la g�n�ration du fichier
        for y=1:2:ceil(Ymax-Ymin)+18
            for z=1:2:ceil(Zmax-Zmin)+18
                
                Lx=[0,1,1,0,0,1,1,0]*2;
                Ly=[0,0,1,1,0,0,1,1]*2;
                Lz=[0,0,0,0,1,1,1,1]*2;
                
                A=zeros(1,24); B=zeros(1,8);
                for i=1:8
                    A(3*i-2:3*i)=[x+Lx(i)-10+Xmin,...
                        y+Ly(i)-10+Ymin,...
                        z+Lz(i)-10+Zmin];
                    
                    B(i)=Mat3D_final(x+Lx(i),y+Ly(i),z+Lz(i));
                end
                AB=[A,B];
                fprintf(fID2,formatSpec,AB);
            end
        end
    end
    fprintf(fID2,'};');
    fclose('all');
end

end

