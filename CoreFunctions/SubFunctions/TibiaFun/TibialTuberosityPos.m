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
%   - plots : a code, '2' to plot local and global results, 
%             '1' to plot only global results, '0' no plot
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
j=0;

% Variable to handle the distance between the end of the diaphysis to the
% start of the epiphysis. Area where the Tibial Tuberosity is thought to be
% located
AltitudesFract = 0.1 : 0.1 : 0.7;

PtsMedThird = zeros(length(AltitudesFract),3);
PtsMid = zeros(length(AltitudesFract),3);
PtsTT = zeros(length(AltitudesFract),3);

% Calculate the position of the tibial tuberosity
for xAlt = AltitudesFract
    j=j+1;
    
    % Find the border of the Proximal Tibia at the altitude xAlt along the
    % Tibia Z axis
    Curve = TriPlanIntersect(ProxTib,CS.Z,-(xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph));
    CSPts0 = Curve.Pts;
    
    % Move the intersection in the Tibia Coordinate Frame
    CSPtsRt = transpose(CS.V'*CSPts0');
    [ Centroid, ~ ] = PlanPolygonCentroid3D( CSPtsRt );
    
    CSPtsRtC0 = bsxfun(@minus,CSPtsRt,Centroid);
    
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
    if fitresultPos.c1 < fitresultNeg.c1
        IDX0 = knnsearch(CSPtsRtC0,[fitresultPos.a1 fitresultPos.b1 0]);
        LegSide = 'R';
    else
        IDX0 = knnsearch(CSPtsRtC0,[-fitresultNeg.a1 fitresultNeg.b1 0]);
        LegSide = 'L';
    end
    
    PtTT0 = CSPtsRtC0(IDX0,:);
    PtTT00 = PtTT0;
    UTT0 =  PtTT0'/norm(PtTT0);
    
    
    
    PtsOnTT0 = unique(CSPtsRtC0(CSPtsRtC0*UTT0>0.925*norm(PtTT0),:),'rows');
    
    PtsOnTTEnd = PtsOnTT0;
    IDX=10^6;
    
    while abs(IDX-IDX0)>1
        IDX0 = IDX;
        
%       [~,IA,~] = intersect(CSPtsRtC0(1:end-1,:),PtsOnTT0,'rows','stable');
        [~,IA,~] = intersect(CSPtsRtC0,PtsOnTT0,'rows','stable');
%       IA = unique(IA,'rows','stable');
        WeightPtsTT = Weights(IA,:);
        
        WeightPtsTT = WeightPtsTT/sum(WeightPtsTT);
        PtMeanTT = sum([WeightPtsTT WeightPtsTT WeightPtsTT].*PtsOnTT0,1);
        
        IDX = knnsearch(CSPtsRtC0,PtMeanTT);
        PtTT = CSPtsRtC0(IDX,:);
        UTT = PtTT'/norm(PtTT);
        PtsOnTT = unique(CSPtsRtC0(CSPtsRtC0*UTT>0.90*norm(PtTT),:),'rows','stable');
        
        PtsOnTTEnd  = unique([PtsOnTTEnd;PtsOnTT],'rows','stable');
        
        PtsOnTT0 = PtsOnTT;
    end
    
    % Finalize the points that are on the TT
    PtsOnTT = PtsOnTTEnd;
    [~,IA,~] = intersect(CSPtsRtC0,PtsOnTT,'rows','stable');
    WeightPtsTT = Weights(IA,:);
    WeightPtsTT = WeightPtsTT/sum(WeightPtsTT);
    PtMeanTT = sum([WeightPtsTT WeightPtsTT WeightPtsTT].*PtsOnTT,1);
    
    IDX = knnsearch(CSPtsRtC0,PtMeanTT);
    PtTT = CSPtsRtC0(IDX,:);
    UTT = PtTT'/norm(PtTT);
    
    % Identify the medial third of the Tibial Tuberosity :
    % Make it medial to lateral like for the Tibia CS
    VTT = [UTT(2) -UTT(1) 0]'; VTT = sign(VTT(2))*VTT;
    [~,IMedialPt] = max(PtsOnTT*VTT);
    [~,ILateralPt] = min(PtsOnTT*VTT);
    
    MedialThirdOfTT = 0.3*PtsOnTT(ILateralPt,:) + 0.7*PtsOnTT(IMedialPt,:);
    MiddleOfTT = 1/2*PtsOnTT(ILateralPt,:) + 1/2*PtsOnTT(IMedialPt,:);
    
    LineEndingOfTT = PtsOnTT(IMedialPt,:)-PtsOnTT(ILateralPt,:);
    VLineEndingOfTT = LineEndingOfTT'/norm(LineEndingOfTT);
    
    % The medial third point
    [~,IDMedThird] = min(abs(bsxfun(@minus,PtsOnTT,MedialThirdOfTT)*VLineEndingOfTT));
    PtMedialThirdOfTT = PtsOnTT(IDMedThird,:);
    
    % The middle of TT
    [~,IDMiddle] = min(abs(bsxfun(@minus,PtsOnTT,MiddleOfTT)*VLineEndingOfTT));
    PtMiddleofTT = PtsOnTT(IDMiddle,:);
    
    if plots >=2
        % Plot fits with data.
        %     figure( 'Name', 'FitParts' );
        %     h = plot( fitresultNeg, xDataNeg, yDataNeg );
        %     hold on
        %     h = plot( fitresultPos, xDataPos, yDataPos , 'm.' );
        %     xlabel X
        %     ylabel Y
        %     grid on
        %     axis equal
        
        
        figure()
        pl3t(CSPtsRtC0,'k-')
        hold on
        pl3t(PtTT00,'g*')
        pl3t(PtTT,'rs')
        pl3t(PtMedialThirdOfTT,'ko')
        pl3t(PtsOnTT,'b.')
        pl3t(MedialThirdOfTT,'ro')
        
        pl3t([0 0 0],'ko')
        
    end
    
    PtsMedThird(j,:) = transpose(CS.V*(PtMedialThirdOfTT + Centroid)') ;
    PtsMid(j,:) = transpose(CS.V*(PtMiddleofTT + Centroid)') ;
    PtsTT(j,:) = transpose(CS.V*(PtTT + Centroid)') ;
    
end

IDPtMedialThird = knnsearch(ProxTib.Points,mean(PtsMedThird));
PtMedialThirdOfTT = ProxTib.Points(IDPtMedialThird,:);

IDPtMiddle = knnsearch(ProxTib.Points,mean(PtsMid));
PtMiddleOfTT = ProxTib.Points(IDPtMiddle,:);

if plots >=1
    figure()
    trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],...
            'FaceAlpha',1,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    axis equal
    light('Position',CS.Origin' + 300*CS.Y + 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 200*CS.Y - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 50*CS.Y + 50*CS.X - 500*CS.Z,'Style','local')
    plotDot( PtsMedThird, 'r', 1.25 )
    plotDot( PtMedialThirdOfTT, 'g', 2.5 )
    plotDot( PtsTT, 'b', 1.25 )
    hold on
    grid off
    lighting gouraud
end

end