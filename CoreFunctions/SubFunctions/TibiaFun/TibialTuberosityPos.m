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

clearvars PtsMedThird circularity Alt
j=0;

% Variable to handle the distance between the end of the diaphysis to the
% start of the epiphysis. Area where the Tibial Tuberosity is thought to be
% located
AltitudesFract = -0.5 : 0.05 : 0.75;


for xAlt = AltitudesFract
    j=j+1;
    Curve = TriPlanIntersect(ProxTib,CS.Z,(xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph));
    CSPts0 = Curve.Pts;
    CSPts0 = interparc(201,CSPts0(:,1),CSPts0(:,2),CSPts0(:,3),'linear');
    CSPtsRt = transpose(CS.V'*CSPts0');
    
    Perimeter = sum(sqrt(sum(diff(CSPtsRt,1).^2,2)));
    Area = polyarea(CSPtsRt(:,1),CSPtsRt(:,2));
    
    circularity(j) = 2 * sqrt(pi) * sqrt(Area)/Perimeter;
    Alt(j) = (xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph);
    
end

[~,I_Cir_Min] = min(circularity);
xAlt0 = AltitudesFract(I_Cir_Min);

AltitudesFract = xAlt0;
j=0;


PtsMedThird = zeros(length(AltitudesFract),3);
PtsMid = zeros(length(AltitudesFract),3);
PtsTT = zeros(length(AltitudesFract),3);

% Calculate the position of the tibial tuberosity
for xAlt = AltitudesFract
    j=j+1;
    
    % Find the border of the Proximal Tibia at the altitude xAlt along the
    % Tibia Z axis
    Curve = TriPlanIntersect(ProxTib,CS.Z,(xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph));
    CSPts0 = Curve.Pts;
    CSPts0 = interparc(201,CSPts0(:,1),CSPts0(:,2),CSPts0(:,3),'linear');
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
    PtTT00 = PtTT0;
    UTT0 =  PtTT0'/norm(PtTT0);
    
    % Renumber curve from the opposite point to the TTA
    [~,Imin] = min(CSPtsRtC0*UTT0);
    CSPtsRtC0 = [CSPtsRtC0(Imin:end-1,:);CSPtsRtC0(1:Imin,:)];
    N = LineNormals2D(CSPtsRtC0(:,1:2));
    
    
    % Get the curvature and keep only the curvature peaks that are on the 
    % anterior face and that make a angle superior to 30° with UTT0
    
    k = LineCurvature2D(CSPtsRtC0(:,1:2));
    
    [pks,locs] = findpeaks(-k);
    
    I1 = find( CSPtsRtC0(locs,:)*UTT0>0.25*max(CSPtsRtC0*UTT0) &... %0.1*max(CSPtsRtC0*UTT0)
               N(locs,:)*UTT0(1:2,:)<sqrt(3)/2 &...
               abs(N(locs,:)*[0;sign(UTT0(1))])>0.5 &...
               pks>quantile(pks,2/3)...
               );
%     %
%                         N(locs,:)*UTT0(1:2,:)<sqrt(3)/2 &...
%               abs(N(locs,:)*[0;sign(UTT0(1))])>0.5 &... %*[0;sign(UTT0(1))])>0.5
%           
    locs_ok = locs(I1);
    
    
    [pks,locs] = findpeaks(k);
    I2 = find( CSPtsRtC0(locs,:)*UTT0>0.33*max(CSPtsRtC0*UTT0) &...
              abs(N(locs,:)*[sign(UTT0(1));0])>0 &...
              pks>quantile(pks,2/3)...
              );
    locs_ok2 = locs(I2);
    
    % Get rid of peaks located outside the TT
    locs_ok2(locs_ok2<locs_ok(1))=[];
    locs_ok2(locs_ok2>locs_ok(end))=[]; 
    
    % Find the Start Antero-Lateral and End Antero-Medial of the TT
    I = find(locs_ok<locs_ok2(1));
    Loc_TT_AL = round((locs_ok(I(end))+locs_ok2(1))/2);
    Loc_TT_AM = min(locs_ok2(end),min(locs_ok(locs_ok>median(locs_ok2))));

    % Find the middle and third of the line separating the TT from the rest
    % of the bone cross section curve
    MedialThirdOfTT = 1/3*CSPtsRtC0(Loc_TT_AL,:) + 2/3*CSPtsRtC0(Loc_TT_AM,:);
    MiddleOfTT = 1/2*CSPtsRtC0(Loc_TT_AL,:) + 1/2*CSPtsRtC0(Loc_TT_AM,:);
    
    % update the Tibial tuberosity direction from 
    UTT0 = MiddleOfTT'/norm(MiddleOfTT);
    
    % Get the separating line direction
    LineEndingOfTT = CSPtsRtC0(Loc_TT_AM,:)-CSPtsRtC0(Loc_TT_AL,:);
    VLineEndingOfTT = LineEndingOfTT'/norm(LineEndingOfTT);
    
    
    % Keep points on the anterior quarter of the cross section curves
    CSPtsRtC0_p = CSPtsRtC0(CSPtsRtC0*UTT0>0.5*max(CSPtsRtC0*UTT0),:);
    
    % The medial third point on curve
    [~,IDMedThird] = min(abs(bsxfun(@minus,CSPtsRtC0_p,MedialThirdOfTT)*VLineEndingOfTT));
    PtMedialThirdOfTT = CSPtsRtC0_p(IDMedThird,:);
    
    % The middle of TT on curve
    [~,IDMiddle] = min(abs(bsxfun(@minus,CSPtsRtC0_p,MiddleOfTT)*VLineEndingOfTT));
    PtMiddleofTT = CSPtsRtC0_p(IDMiddle,:);
        
    if plots ==1
        plotCurvature2D(CSPtsRtC0,k)
        hold on
        pl3t(CSPtsRtC0(1,:),'ko')
        pl3t(PtTT0,'ko')
        pl3t(CSPtsRtC0(locs_ok2,:)','ro')
        pl3t(CSPtsRtC0(locs_ok,:)','ks')
        pl3t(mean(CSPtsRtC0(locs_ok2(1):locs_ok2(end),:)),'m+')
        pl3t(CSPtsRtC0(Loc_TT_AL,:),'mo')
        pl3t([CSPtsRtC0(Loc_TT_AL,:);CSPtsRtC0(Loc_TT_AL,:);...
            CSPtsRtC0(Loc_TT_AM,:);CSPtsRtC0(Loc_TT_AM,:)],'k-')
        pl3t(PtMiddleofTT,'k*')
        pl3t(PtMedialThirdOfTT,'r*')
        pl3t(zeros(3,5),'ko')
    end
    
    PtsMedThird(j,:) = transpose(CS.V*(PtMedialThirdOfTT + Centroid)') ;
    PtsMid(j,:) = transpose(CS.V*(PtMiddleofTT + Centroid)') ;
%     PtsTT(j,:) = transpose(CS.V*(PtTT + Centroid)') ;
    
Perimeter = sum(sqrt(sum(diff(CSPtsRtC0,1).^2,2)));
Area = polyarea(CSPtsRtC0(:,1),CSPtsRtC0(:,2));

circularity(j) = 2 * sqrt(pi) * sqrt(Area)/Perimeter;
Alt(j) = (xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph);
end

% IDPtMedialThird = knnsearch(ProxTib.Points,mean(PtsMedThird));
% PtMedialThirdOfTT = ProxTib.Points(IDPtMedialThird,:);
% 
% IDPtMiddle = knnsearch(ProxTib.Points,mean(PtsMid));
% PtMiddleOfTT = ProxTib.Points(IDPtMiddle,:);


IDPtMedialThird = knnsearch(ProxTib.Points,PtsMedThird);
PtMedialThirdOfTT = ProxTib.Points(IDPtMedialThird,:);

IDPtMiddle = knnsearch(ProxTib.Points,PtsMid);
PtMiddleOfTT = ProxTib.Points(IDPtMiddle,:);

if plots ==1
    figure()
    trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],...
            'FaceAlpha',1,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    axis equal
    light('Position',CS.Origin' + 500*CS.Y + 500*CS.X - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 500*CS.Y - 500*CS.X,'Style','local')
    light('Position',CS.Origin' - 500*CS.Y + 500*CS.X + 500*CS.Z,'Style','local')
    plotDot( PtsMedThird, 'r', 1.25 )
    plotDot( PtMedialThirdOfTT, 'g', 2.5 )
    plotDot( PtsMid, 'b', 1.25 )
%     plotDot( PtsTT, 'b', 1.25 )
    hold on
    grid off
    lighting gouraud
    
    figure()
    plot(Alt,circularity,'.-')
end

end