function [ Results ] = TibiaCS( ProxTib , DistTib)
% Fit an ACS on a Tibia composed of the proximal tibia and the tibial part of
% the ankle

addpath(strcat(pwd,'\SubFunctions\SurfFit'));
addpath(strcat(pwd,'\SubFunctions\TibiaFun'));
addpath(strcat(pwd,'\SubFunctions\TriFun'));

% Unite both distal and proximal tibia mesh
Tibia = TriUnite(ProxTib,DistTib);

[ InertiaMatrix, Center ] = TriInertiaProperties( Tibia );
[V_all,~] = eig(InertiaMatrix);
Center0 = Center;


% Initial estimate of the Inf-Sup axis Z0 - Check that the distal tibia
% is 'below': the proximal tibia, invert Z0 direction otherwise;
Z0 = V_all(:,1);
Z0 = sign((mean(ProxTib.Points)-mean(DistTib.Points))*Z0)*Z0;

Minertia = V_all;

%% Distal Tibia

% 1st : Find triangles with less than 30° relative to the tibia principal
% inertia axis (longitudinal) and low curvature then smooth the result with close
% morphology operation

Alt =  min(DistTib.Points*Z0)+1 : 0.3 : max(DistTib.Points*Z0)-1;
Area=[];
for d = Alt
    [ ~ , Area(end+1), ~ ] = TriPlanIntersect( DistTib, Z0 , d );
end
[~,Imax] = max(Area);
[ Curves ,  ~ , ~ ] = TriPlanIntersect( DistTib, Z0 , Alt(Imax) );

CenterAnkleInside = PlanPolygonCentroid3D( Curves.Pts);


[Cmean,~,~,~,~,~]=TriCurvature(DistTib,false);

AnkleArtSurfNodesOK0 =  find(Cmean>quantile(Cmean,0.65) & ...
    Cmean<quantile(Cmean,0.95) & ...
    rad2deg(acos(-DistTib.vertexNormal*Z0))<30);

AnkleArtSurf0 = TriReduceMesh(DistTib,[],double(AnkleArtSurfNodesOK0));
AnkleArtSurf0 = TriCloseMesh(DistTib,AnkleArtSurf0,6);
AnkleArtSurf0 = TriOpenMesh(DistTib,AnkleArtSurf0,4);
AnkleArtSurf0 = TriConnectedPatch( AnkleArtSurf0, mean(AnkleArtSurf0.Points));

% 2nd : fit a polynomial surface to it AND Exclude points that are two far (1mm) from the fitted surface,
% then smooth the results with open & close morphology operations
TibiaElmtsIDOK = AnkleSurfFit( AnkleArtSurf0, DistTib, V_all );
AnkleArtSurf = TriReduceMesh(DistTib , TibiaElmtsIDOK);
AnkleArtSurf = TriErodeMesh(AnkleArtSurf,2);
AnkleArtSurf = TriCloseMesh(DistTib,AnkleArtSurf,6);
AnkleArtSurf = TriOpenMesh(DistTib,AnkleArtSurf,4);
AnkleArtSurf = TriCloseMesh(DistTib,AnkleArtSurf,2);


% Filter Elmts that are not
AnkleArtSurfProperties = TriMesh2DProperties( AnkleArtSurf );
ZAnkleSurf = AnkleArtSurfProperties.meanNormal;
CAnkle = AnkleArtSurfProperties.onMeshCenter;
AnkleArtSurfElmtsOK = find(rad2deg(acos(AnkleArtSurf.faceNormal*ZAnkleSurf))<35 & ...
    sum(bsxfun(@minus,AnkleArtSurf.incenter,CAnkle).*AnkleArtSurf.faceNormal,2)./...
    sqrt(sum(bsxfun(@minus,AnkleArtSurf.incenter,CAnkle).^2,2))<0.1);
AnkleArtSurf = TriReduceMesh(AnkleArtSurf,AnkleArtSurfElmtsOK);
AnkleArtSurf = TriOpenMesh(DistTib,AnkleArtSurf,4);
AnkleArtSurf = TriCloseMesh(DistTib,AnkleArtSurf,10);


% Method : 1) Fit a LS plan to the art surface, 2) Inset the plan 5mm
% 3) Get the center of Shape with shape intersection 4) Project Center Back
% to original plan


[nAAS,dAAS] = LS_Plan(AnkleArtSurf.Points);
nAAS = sign(nAAS'*Z0)*nAAS; dAAS = sign(nAAS'*Z0)*dAAS;

[ Curves , ~, ~ ] = TriPlanIntersect( DistTib, nAAS , -dAAS+5 );

Centr = PlanPolygonCentroid3D( Curves.Pts );

CenterAnkle = Centr -5*nAAS';
%% Find a pseudo medioLateral Axis :
% Most Distal point of the medial malleolus (MDMMPt)
ZAnkleSurf = AnkleArtSurfProperties.meanNormal;
[~,I] = max(DistTib.Points*ZAnkleSurf);
MDMMPt = DistTib.Points(I,:);

U_tmp = MDMMPt - CenterAnkle;
Y0 = U_tmp' - (U_tmp*Z0)*Z0; Y0=Y0/norm(Y0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  END OF DISTAL TIBIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Proximal Tibia

% First 0.5 mm in Start and End are not accounted for, for stability.

Alt = linspace( min(ProxTib.Points*Z0)+0.5 ,max(ProxTib.Points*Z0)-0.5, 100);
Area = zeros(size(Alt));
i=0;
for d = Alt
    i=i+1;
    [ ~ , Area(i), ~ ] = TriPlanIntersect( ProxTib, Z0 , d );
end
AltAtMax = Alt(Area==max(Area));

[~,Zepi,~] = EpiphysisDiaphysisEnds(Alt, Area);

ElmtsEpi = find(ProxTib.incenter*Z0>Zepi); % & rad2deg(acos(ProxTib.faceNormal*Z0))<45;
EpiTib = TriReduceMesh( ProxTib, ElmtsEpi );


[Cmean,Cgaussian,~,~,~,~]=TriCurvature(EpiTib,false);

Curvtr = sqrt(4*Cmean.^2-2*Cgaussian);

NodesEpiArtSurfOK = find(rad2deg(acos(EpiTib.vertexNormal*Z0))<35 &...
    Curvtr<quantile(Curvtr,0.25)) ;
Pcondyle = EpiTib.Points(NodesEpiArtSurfOK,:);
EpiTibArt = TriReduceMesh( EpiTib, [] , NodesEpiArtSurfOK );
EpiTibArt = TriCloseMesh( EpiTib, EpiTibArt, 6 );
[Ztp,d] = LS_Plan(Pcondyle); Ztp = sign(Z0'*Ztp)*Ztp;

% Associate an ellipse to the kept elements external border
[ Xel, Yel, ellipsePts , ellipsePpties] = EllipseOnEdge( EpiTibArt, Ztp , d );
a = ellipsePpties.a;
b = ellipsePpties.b;

% Define 3 points on the Medial Articular Surface and 3 on the Lateral
MedPtsInit = mean(ellipsePts) + 2/3*b*Yel';
MedPtsInit = [MedPtsInit; MedPtsInit - 1/3*a*Xel'; MedPtsInit + 1/3*a*Xel'];
LatPtsInit = mean(ellipsePts) - 2/3*b*Yel';
LatPtsInit = [LatPtsInit; LatPtsInit - 1/3*a*Xel'; LatPtsInit + 1/3*a*Xel'];

% Use those points to isolate a patch of elements that are on the  medial
% and lateral articular surface
EpiTibArtMed = TriConnectedPatch( EpiTibArt, MedPtsInit );
EpiTibArtLat = TriConnectedPatch( EpiTibArt, LatPtsInit );
EpiTibArt = TriUnite(EpiTibArtMed,EpiTibArtLat);


% Repeat the procedure to update the ellipse
[Ztp,d] = LS_Plan(EpiTibArt.Points); Ztp = sign(Z0'*Ztp)*Ztp;
[ Xel, Yel, ellipsePts , ellipsePpties] = EllipseOnEdge( EpiTibArt, Ztp , d );
a = ellipsePpties.a;
b = ellipsePpties.b;

% Defined an hourglass shaped elements set of the tibial plateau as not 
% containing articular surfaces elements : EpiTibCenterRidge
Xel = sign(Xel'*Y0)*Xel;
Yel = sign(Yel'*Y0)*Yel;

d = mean(ellipsePts)*Xel+0.5*a;

[ Curves , ~, ~ ] = TriPlanIntersect( ProxTib, Xel , d );

MedialPts_tmp = Curves(1).Pts(bsxfun(@minus,Curves(1).Pts,mean(ellipsePts))*Yel>0,:);
[~,IDPtsMax] = max(MedialPts_tmp*Z0);
PtsMax = MedialPts_tmp(IDPtsMax,:); %+mean(ellipsePts);

U_tmp =  transpose(PtsMax-mean(ellipsePts));

np = cross(U_tmp,Ztp); np = sign(cross(Xel,Yel)'*Z0)*np/norm(np);
dp = mean(ellipsePts)*np;

nm = Yel;
dm = mean(ellipsePts)*nm;

NodesOnCenterID = find(sign(EpiTib.Points*np-dp) + sign(EpiTib.Points*nm-dm)>0.1);
EpiTibCenterRidgeMed = TriReduceMesh( EpiTib, [] , NodesOnCenterID );


LateralPts_tmp = Curves(1).Pts(bsxfun(@minus,Curves(1).Pts,mean(ellipsePts))*Yel<0 & ...
    bsxfun(@minus,Curves(1).Pts,mean(ellipsePts))*Yel>-b/3&...
    abs(bsxfun(@minus,Curves(1).Pts,mean(ellipsePts))*Z0)<a/2,:);
[~,IDPtsMax] = min(LateralPts_tmp*Z0);
PtsMax = LateralPts_tmp(IDPtsMax,:); %+mean(ellipsePts);

U_tmp =  transpose(PtsMax-mean(ellipsePts));

np = cross(U_tmp,Ztp); np = -sign(cross(Xel,Yel)'*Z0)*np/norm(np);
dp = mean(ellipsePts)*np;

nm = -Yel;
dm = mean(ellipsePts)*nm;

NodesOnCenterID = find(sign(EpiTib.Points*np-dp) + sign(EpiTib.Points*nm-dm)>0.1);
EpiTibCenterRidgeLat = TriReduceMesh( EpiTib, [] , NodesOnCenterID );
EpiTibCenterRidgeLat = TriDilateMesh(EpiTib, EpiTibCenterRidgeLat,5);


EpiTibCenterRidge = TriUnite(EpiTibCenterRidgeLat,EpiTibCenterRidgeMed);

% Update the previously identified AS elements to remove the one 
% (potentially) located on the central part of the tibial plateau
MedPtsInit = mean(ellipsePts) + 2/3*b*Yel';
MedPtsInit = [MedPtsInit; MedPtsInit - 1/3*a*Xel'; MedPtsInit + 1/3*a*Xel'];
LatPtsInit = mean(ellipsePts) - 2/3*b*Yel';
LatPtsInit = [LatPtsInit; LatPtsInit - 1/3*a*Xel'; LatPtsInit + 1/3*a*Xel'];


EpiTibArt = TriDifferenceMesh(EpiTibArt,EpiTibCenterRidge);
EpiTibArt = TriDifferenceMesh(EpiTibArt,EpiTibCenterRidge);
EpiTibArtMed = TriConnectedPatch( EpiTibArt, MedPtsInit );
EpiTibArtLat = TriConnectedPatch( EpiTibArt, LatPtsInit );
EpiTibArt = TriUnite(EpiTibArtMed,EpiTibArtLat);
EpiTibArt = TriOpenMesh(EpiTib,EpiTibArt, 15);
EpiTibArt = TriCloseMesh(EpiTib,EpiTibArt, 30);

% Update the plan and the ellipse
[Ztp,d] = LS_Plan(EpiTibArt.Points); Ztp = sign(Z0'*Ztp)*Ztp;
[ Xel, Yel, ellipsePts , ellipsePpties] = EllipseOnEdge( EpiTibArt, Ztp , d );
a = ellipsePpties.a;
b = ellipsePpties.b;
Xel = sign(Xel'*Y0)*Xel;
Yel = sign(Yel'*Y0)*Yel;
MedPtsInit = mean(ellipsePts) + 2/3*b*Yel';
MedPtsInit = [MedPtsInit; MedPtsInit - 1/3*a*Xel'; MedPtsInit + 1/3*a*Xel'];
LatPtsInit = mean(ellipsePts) - 2/3*b*Yel';
LatPtsInit = [LatPtsInit; LatPtsInit - 1/3*a*Xel'; LatPtsInit + 1/3*a*Xel'];


EpiTibArtMed = TriConnectedPatch( EpiTibArt, MedPtsInit);
EpiTibArtLat = TriConnectedPatch( EpiTibArt, LatPtsInit );

EpiTibArtMedElmtsOK = find(abs(EpiTibArtMed.incenter*Ztp+d)<5 & ...
    EpiTibArtMed.faceNormal*Ztp>0.9 );
EpiTibArtMed = TriReduceMesh(EpiTibArtMed,EpiTibArtMedElmtsOK);
EpiTibArtMed = TriOpenMesh(EpiTib,EpiTibArtMed,2);
EpiTibArtMed = TriConnectedPatch( EpiTibArtMed, MedPtsInit );
EpiTibArtMed = TriCloseMesh(EpiTib,EpiTibArtMed,10);

EpiTibArtLatElmtsOK = find(abs(EpiTibArtLat.incenter*Ztp+d)<5 & ...
    EpiTibArtLat.faceNormal*Ztp>0.9 );
EpiTibArtLat = TriReduceMesh(EpiTibArtLat,EpiTibArtLatElmtsOK);
EpiTibArtLat = TriOpenMesh(EpiTib,EpiTibArtLat,2);
EpiTibArtLat = TriConnectedPatch( EpiTibArtLat, LatPtsInit );
EpiTibArtLat = TriCloseMesh(EpiTib,EpiTibArtLat,10);

EpiTibArt = TriUnite(EpiTibArtMed,EpiTibArtLat);

[Ztp,d] = LS_Plan(EpiTibArt.Points); Ztp = sign(Z0'*Ztp)*Ztp;
[ Xel, Yel, ~ , ellipsePpties] = EllipseOnEdge( EpiTibArt, Ztp , d );
a = ellipsePpties.a;
b = ellipsePpties.b;
Xel = sign(Xel'*Y0)*Xel;
Yel = sign(Yel'*Y0)*Yel;

EpiTibArtMedElmtsOK = find(abs(EpiTibArtMed.incenter*Ztp+d)<5 & ...
    EpiTibArtMed.faceNormal*Ztp>0.95 );
EpiTibArtMed = TriReduceMesh(EpiTibArtMed,EpiTibArtMedElmtsOK);
EpiTibArtMed = TriOpenMesh(EpiTib,EpiTibArtMed,2);
EpiTibArtMed = TriConnectedPatch( EpiTibArtMed, MedPtsInit );
EpiTibArtMed = TriCloseMesh(EpiTib,EpiTibArtMed,10);

EpiTibArtLatElmtsOK = find(abs(EpiTibArtLat.incenter*Ztp+d)<3 & ...
    EpiTibArtLat.faceNormal*Ztp>0.95 );
EpiTibArtLat = TriReduceMesh(EpiTibArtLat,EpiTibArtLatElmtsOK);
EpiTibArtLat = TriOpenMesh(EpiTib,EpiTibArtLat,2);
EpiTibArtLat = TriConnectedPatch( EpiTibArtLat, LatPtsInit );
EpiTibArtLat = TriCloseMesh(EpiTib,EpiTibArtLat,10);

EpiTibArt = TriUnite(EpiTibArtMed,EpiTibArtLat);

[Ztp,d] = LS_Plan(EpiTibArt.Points); Ztp = sign(Z0'*Ztp)*Ztp; d = sign(Z0'*Ztp)*d;
[ Xel, Yel, ~ , ~] = EllipseOnEdge( EpiTibArt, Ztp , d );
Xel = sign(Xel'*Y0)*Xel;
Yel = sign(Yel'*Y0)*Yel;


%% Technic PIAASL : Compute the inertial axis of a slice of the tibial plateau
% 10% below and the 5% above : Fill it with equally spaced points to
% simulate inside volume
%

H = 0.1 * sqrt(4*0.75*max(Area)/pi); % 10% of the general size of the TP

Alt_TP = linspace( -d-H ,-d+0.5*H, 20);
PointSpace = mean(diff(Alt_TP));
TPLayerPts = zeros(round(length(Alt_TP)*1.1*max(Area)/PointSpace^2),3);
j=0;
for alt = Alt_TP
    [ Curves , ~ , ~ ] = TriPlanIntersect( EpiTib, Ztp , alt );
    for c=1:length(Curves)
        
        Pts_Tmp = Curves(c).Pts*[Xel Yel Ztp];
        xmg = min(Pts_Tmp(:,1)) -0.1 : PointSpace : max(Pts_Tmp(:,1)) +0.1 ;
        ymg = min(Pts_Tmp(:,2)) -0.1 : PointSpace : max(Pts_Tmp(:,2)) +0.1;
        [XXmg , YYmg] = meshgrid(xmg,ymg);
        in = inpolygon(XXmg(:),YYmg(:),Pts_Tmp(:,1),Pts_Tmp(:,2));
        Iin = find(in, 1);
        if ~isempty(Iin)
            i = j+1;
            j=i+length(find(in))-1;
            TPLayerPts(i:j,:) = transpose([Xel Yel Ztp]*[XXmg(in),YYmg(in),ones(length(find(in)),1)*alt]');
        end
    end
    
end

TPLayerPts(j+1:end,:) = [];

[V,~] = eig(cov(TPLayerPts));

Xtp = V(:,2); Ytp = V(:,3);
Xtp = sign(Xtp'*Y0)*Xtp;
Ytp = sign(Ytp'*Y0)*Ytp;

idx = kmeans(TPLayerPts,2);

[ CenterMed ] = ProjectOnPlan( mean(TPLayerPts(idx==1,:)) , Ztp , d );
[ CenterLat ] = ProjectOnPlan( mean(TPLayerPts(idx==2,:)) , Ztp , d );

CenterKnee = 0.5*( CenterMed + CenterLat);

Zmech = CenterKnee - CenterAnkleInside; Zmech = Zmech' / norm(Zmech);

% Final ACS
Xend = cross(Ytp,Zmech)/norm(cross(Ytp,Zmech));
Yend = cross(Zmech,Xend);

Yend = sign(Yend'*Y0)*Yend;
Zend = Zmech;
Xend = cross(Yend,Zend);

Vend = [Xend Yend Zend];

% Result write
Results.CenterAnkle = CenterAnkle;
Results.CenterKnee = CenterKnee;
Results.Ztp = Ztp;
Results.Ytp = Ytp;
Results.Xtp = Xtp;
Results.PlanTPd = d;

Results.Xend = Xend;
Results.Yend = Yend;
Results.Zend = Zend;
Results.Vend = Vend;

%% Inertia Results
Yi = V_all(:,2); Yi = sign(Yi'*Y0)*Yi;
Xi = cross(Yi,Z0);


Results.Center0 = Center0;
Results.Zinertia = Z0;
Results.Yinertia = Yi;
Results.Xinertia = Xi;
Results.Minertia = [Xi Yi Z0];


end

