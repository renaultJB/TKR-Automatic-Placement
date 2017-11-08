function [ PtMedialThirdOfTT, LegSide, PtsMedThird, PtsTT, PtMiddleOfTT ] = TibialTuberosityPos(ProxTib, CS , plots)
%--------------------------------------------------------------------------
% TibialTuberosityPos : Find a Point on medial third of the tibial
% tubercule
% Paired with CSTibia.m the script can identify id the tibia is from a
% right or a left leg
%
%   Lexic :
%       TT ? Tubial Tubercule
%       Med ? Medial
%       CS ? Coordinate system
%   Inputs:
%   - Proxtib : Triangulation object representing the tibia bone
%   - CS : A structure of the constructed coordinate system associated to 
%          the tibia bone
%   - plots : a boolean, '1' to plot results, '0' to not plot
%
%   Outputs:
%   - PtMedialThirdOfTT : A point on the original mesh that is between the
%                         medial and center thirds of the tibial tuberosity
%   - Legside : string that describe the leg side either 'L' or 'R'
%   - PtsMedThird : A Nx3 matrix of the points use to get PtMedialThirdOfTT
%   - PtsTT : A Nx3 matrix of the points located on the middle of the TTA
%   - PtMiddleOfTT : Same as PtMedialThirdOfTT but located on the middle
%                    of the TTA
%
%--------------------------------------------------------------------------
if nargin < 3
    plots = 0;
end

clearvars PtsMedThird


% Variable to handle the distance between the end of the diaphysis to the
% start of the epiphysis. Area where the Tibial Tuberosity is thought to be
% located

AltitudesFract = 0;

PtsMedThird = zeros(length(AltitudesFract),3);
PtsMid = zeros(length(AltitudesFract),3);
PtsTT = zeros(length(AltitudesFract),3);

% Calculate the position of the tibial tuberosity
for xAlt = AltitudesFract
    
    % Find the border of the Proximal Tibia at the altitude xAlt along the
    % Tibia Z axis
    Curve = TriPlanIntersect(ProxTib,CS.Z,(xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph));
    CSPts0 = Curve.Pts;
    Perimeter = sum(sqrt(sum(diff(CSPts0,1).^2,2)));
    CSPts0 = interparc(round(Perimeter),CSPts0(:,1),CSPts0(:,2),CSPts0(:,3),'linear');
    CSPtsRt = transpose(CS.V'*CSPts0');
    [ Centroid, ~ ] = PlanPolygonCentroid3D( CSPtsRt );
    CSPtsRtC0 = bsxfun(@minus,CSPtsRt,Centroid);
    
    % Check that the normals are pointing outside, otherwise invert
    % ordering of vertex
    N = LineNormals2D(CSPtsRtC0(:,1:2));
    OffSetCS = [CSPtsRtC0(:,1)+5*N(:,2),CSPtsRtC0(:,2)+5*N(:,2)];
    InSetCS = [CSPtsRtC0(:,1)-5*N(:,2),CSPtsRtC0(:,2)-5*N(:,2)];
    
    if polyarea(OffSetCS(:,1),OffSetCS(:,2)) < polyarea(InSetCS(:,1),InSetCS(:,2))
        CSPtsRtC0 = CSPtsRtC0(end:-1:1,:);
        CSPtsRt = CSPtsRt(end:-1:1,:);
        N = LineNormals2D(CSPtsRtC0(:,1:2));
    end
        
    %Calculate the weight of each vertex of the curves,
    % the weight is the sum of the length two connected edges
    
    Lim1 = sqrt(sum(diff([CSPtsRtC0(end-1,:) ; CSPtsRtC0(1:end-1,:)],1,1).^2,2));
    Li = sqrt(sum(diff(CSPtsRtC0,1,1).^2,2));
    
    Weights = Lim1 + Li;
    
    Weights(end+1) = Weights(1);
    
    CSPtsPos = CSPtsRtC0(CSPtsRtC0(:,1)>0,:);
    CSPtsNeg = CSPtsRtC0(CSPtsRtC0(:,1)<0,:);
    
    Ypos = CSPtsPos(:,1); Xpos = CSPtsPos(:,2);
    Yneg = -CSPtsNeg(:,1); Xneg = CSPtsNeg(:,2);
    
    
    %% Fit: 'FitNegPart' & 'FitPosPart'.
    [xDataNeg, yDataNeg] = prepareCurveData( Xneg, Yneg );
    [xDataPos, yDataPos] = prepareCurveData( Xpos, Ypos );
    
    % Set up fittype and options.
    ft = fittype( 'gauss1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0];
    opts.StartPoint = [20 0 10];
    
    
    % Fit model to data.
    [fitresultNeg, ~] = fit( xDataNeg, yDataNeg, ft, opts );
    [fitresultPos, ~] = fit( xDataPos, yDataPos, ft, opts );
    
    % Test the most the least flat side(higher sigma value of gaussian fit)
    % and find the closest point on curve to the gaussian fit peak
    if fitresultPos.c1 < fitresultNeg.c1
        IDX0 = knnsearch(CSPtsRtC0,[fitresultPos.a1 fitresultPos.b1 0]);
        LegSide = 'R' ;
    else
        IDX0 = knnsearch(CSPtsRtC0,[-fitresultNeg.a1 fitresultNeg.b1 0]);
        LegSide = 'L' ;
    end
    
    PtTT0 = CSPtsRtC0(IDX0,:);
    UTT0 =  PtTT0'/norm(PtTT0);
    
    %% Little trick so that the code works on both leg side
    if LegSide == 'L'
        CSPtsRtC0(:,1) = -CSPtsRtC0(:,1);
        CSPtsRtC0 = CSPtsRtC0(end:-1:1,:);
        UTT0(1) = -UTT0(1);
        PtTT0(1) = -PtTT0(1);
    end
        
    
    % Renumber curve from the opposite point to the TTA
    [~,Imin] = min(CSPtsRtC0*UTT0);
    CSPtsRtC0 = [CSPtsRtC0(Imin:end-1,:);CSPtsRtC0(1:Imin,:)];
    N = LineNormals2D(CSPtsRtC0(:,1:2));
    
    % Get the curvature and keep only the curvature peaks that are on the 
    % anterior face and that make a angle superior to 30° with UTT0
    
    k = LineCurvature2D(CSPtsRtC0(:,1:2));
    
    [pks,locs] = findpeaks(-k);
    
    % Find the maximal curvature in the Medial and Lateral side of the TTA
    
    ILat = find( CSPtsRtC0(locs,:)*UTT0>0.25*max(CSPtsRtC0*UTT0) &... %0.1*max(CSPtsRtC0*UTT0)
        CSPtsRtC0(locs,:)*[cos(pi/3);sin(pi/3);0]<0);
    pksLat = pks(ILat);
    locsLat= locs(ILat);
    locsLatMax1 = locsLat(pksLat==max(pksLat));
    
    IMed = find( CSPtsRtC0(locs,:)*UTT0>0.25*max(CSPtsRtC0*UTT0) &... %0.1*max(CSPtsRtC0*UTT0)
                 N(locs,2)>0 &...
                 CSPtsRtC0(locs,:)*[cos(pi/3);sin(pi/3);0]>0);
    pksMed = pks(IMed);
    locsMed= locs(IMed);
    locsMedMax1 = locsMed(pksMed==max(pksMed));
    
    
    %    ILat = find( CSPtsRtC0(locs,:)*UTT0>0.25*max(CSPtsRtC0*UTT0) &... %0.1*max(CSPtsRtC0*UTT0)
    [pks,locs] = findpeaks(k);
    
    ILat = find( CSPtsRtC0(locs,:)*UTT0>0.25*max(CSPtsRtC0*UTT0) &... %0.1*max(CSPtsRtC0*UTT0)
                 CSPtsRtC0(locs,:)*[cos(pi/3);sin(pi/3);0]<0);
    pksLat = pks(ILat);
    locsLat= locs(ILat);
    locsLatMax2 = locsLat(pksLat==max(pksLat));
    
    IMed = find( CSPtsRtC0(locs,:)*UTT0>0.33*max(CSPtsRtC0*UTT0) &... %0.1*max(CSPtsRtC0*UTT0)
                 N(locs,2)>0 &...
                 CSPtsRtC0(locs,:)*[cos(pi/3);sin(pi/3);0]>0);
    pksMed = pks(IMed);
    locsMed= locs(IMed);
    locsMedMax2 = locsMed(pksMed==max(pksMed));
    
    
    %% Get the direction of the antero-medial face  
    [~,Imax] = max(CSPtsRtC0(:,2)); % the most medial vertex in CS
    Pts_AM_Face = CSPtsRtC0(locsMedMax2:Imax,:);
    [V_AM_Face,~] = eig(cov(Pts_AM_Face));
    U_AM_Face = V_AM_Face(:,3);
    U_AM_Face = sign(UTT0'*U_AM_Face)*U_AM_Face;
    N_AM_Face = V_AM_Face(:,2);
    N_AM_Face = sign(N_AM_Face(2))*N_AM_Face;
    
    % find the limit of the TTA on the antero-medial face  
    IMed = find(CSPtsRtC0*[cos(pi/3);sin(pi/3);0]>...
                 0.2*max(CSPtsRtC0*[cos(pi/3);sin(pi/3);0]));     
    CSPtsRtC0_M = CSPtsRtC0(IMed,:);
    [~,Loc_TT_AM] = min(abs(bsxfun(@minus,CSPtsRtC0_M,CSPtsRtC0(locsLatMax1,:))*U_AM_Face));
    Loc_TT_AM = IMed(1) + Loc_TT_AM + round(0.666*(locsLatMax2-locsLatMax1));
    Loc_TT_AM = min(Loc_TT_AM,locsMedMax1);
    
    % find the limit of the TTA on the antero-lateral face  
    Loc_TT_AL = round((1*locsLatMax2+3*locsLatMax1)/4);
    
    % Get the TT center and medial third limits
    Loc_TT_Center = round(0.5*Loc_TT_AL + 0.5*Loc_TT_AM);
    Loc_TT_M_Third = round(1/3*Loc_TT_AL + 2/3*Loc_TT_AM);
    
    PtMedialThirdOfTT = CSPtsRtC0(Loc_TT_M_Third,:);
    PtMiddleofTT = CSPtsRtC0(Loc_TT_Center,:);
    
    
    if plots ==1
        plotCurvature2D(CSPtsRtC0,k)
        hold on
        pl3t(zeros(3,5),'ko')
        pl3t(CSPtsRtC0(locsLatMax1,:)','ro')
        pl3t(CSPtsRtC0(locsMedMax1,:)','rs')
        pl3t(CSPtsRtC0(locsLatMax2,:)','bo')
        pl3t(CSPtsRtC0(locsMedMax2,:)','bs')
        pl3t(CSPtsRtC0(Loc_TT_AM,:)','r*')
        pl3t(CSPtsRtC0(Loc_TT_AL,:),'mo')
        pl3t([CSPtsRtC0(Loc_TT_AL,:);CSPtsRtC0(Loc_TT_AL,:);...
            CSPtsRtC0(Loc_TT_AL,:);CSPtsRtC0(Loc_TT_AM,:)],'k-')
        pl3t(PtMiddleofTT,'m*')
        pl3t(PtMiddleofTT,'ms')
        pl3t(PtMedialThirdOfTT,'r*')
        pl3t(PtMedialThirdOfTT,'rs')
    end
    
    if LegSide == 'L'
        PtMedialThirdOfTT(:,1) = -PtMedialThirdOfTT(:,1);
        PtMiddleofTT(1) = -PtMiddleofTT(1);
    end
    
    PtsMedThird = transpose(CS.V*(PtMedialThirdOfTT + Centroid)') ;
    PtsMid = transpose(CS.V*(PtMiddleofTT + Centroid)') ;
    

    

end
IDPtMedialThird = knnsearch(ProxTib.Points,PtsMedThird);
PtMedialThirdOfTT = ProxTib.Points(IDPtMedialThird,:);

IDPtMiddle = knnsearch(ProxTib.Points,PtsMid);
PtMiddleOfTT = ProxTib.Points(IDPtMiddle,:);

if plots ==1
    figure()
    trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],...
            'FaceAlpha',1,'edgecolor','none');
    hold on
    axis equal
    light('Position',CS.Origin' + 500*CS.Y + 500*CS.X - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 500*CS.Y - 500*CS.X,'Style','local')
    light('Position',CS.Origin' - 500*CS.Y + 500*CS.X + 500*CS.Z,'Style','local')
    plotDot( PtsMedThird, 'r', 1.25 )
    plotDot( PtsMid, 'b', 1.25 )
    hold on
    grid off
    lighting gouraud
end

end