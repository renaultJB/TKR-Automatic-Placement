function [ PtMedialThirdOfTT, LegSide, PtsMedThird, PtsTT ] = TibialTuberosityPos(ProxTib, CS , plots)
% TibialTuberosityPos : Find a Point on medial third of the tibial
% tuberosity
% Paired with CSTibia.m the script can identify id the tibia is from a
% right or a left leg
% Lexic :
% TT = Tubial Tuberosity

if nargin < 3
    plots = 0;
end

clearvars PtsMedThird
j=0;

PtsMedThird = zeros(7,3);
PtsTT = zeros(7,3);

% Calculate the position of the
for xAlt = 0.2 : 0.1 : 0.8
    j=j+1;
    
    % Find the border of the Proximal Tibia at the altitude xAlt along the
    % Tibia Z axis
    Curve = TriPlanIntersect(ProxTib,CS.Z,(xAlt*CS.AltStartEpiph+(1-xAlt)*CS.AltEndDiaph));
    CSPts0 = Curve.Pts;
    CSPtsRt = transpose(CS.V'*CSPts0');
    [ Centroid, ~ ] = PlanPolygonCentroid3D( CSPtsRt );
    
    CSPtsRtC0 = bsxfun(@minus,CSPtsRt,Centroid);
    
    %Calculate the weight of each vertex of the curves,
    % the weight is the sum of the length two connected edges
    
    Lim1 = sqrt(sum(diff([CSPtsRtC0(end-1,:) ; CSPtsRtC0(1:end-1,:)],1,1).^2,2));
    Li = sqrt(sum(diff(CSPtsRtC0,1,1).^2,2));
    
    Weights = Lim1 + Li;
    
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
    
    % Test the most the least flat side
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
    
    
    PtsOnTT0 = unique(CSPtsRtC0(CSPtsRtC0*UTT0>0.925*norm(PtTT0),:),'rows');
    
    PtsOnTTEnd = PtsOnTT0;
    IDX=10^6;
    
    while abs(IDX-IDX0)>1
        IDX0 = IDX;
        
        [~,IA,~] = intersect(CSPtsRtC0(1:end-1,:),PtsOnTT0,'rows','stable');
        WeightPtsTT = Weights(IA,:);      
        
        WeightPtsTT = WeightPtsTT/sum(WeightPtsTT);
        PtMeanTT = sum([WeightPtsTT WeightPtsTT WeightPtsTT].*PtsOnTT0,1);
        
        IDX = knnsearch(CSPtsRtC0,PtMeanTT);
        PtTT = CSPtsRtC0(IDX,:);
        UTT = PtTT'/norm(PtTT);
        PtsOnTT = unique(CSPtsRtC0(CSPtsRtC0*UTT>0.90*norm(PtTT),:),'rows');
        
        PtsOnTTEnd  = unique([PtsOnTTEnd;PtsOnTT],'rows','stable');
        
        PtsOnTT0 = PtsOnTT;
    end
    
    % Finalise the points that are on the TT
    PtsOnTT = PtsOnTTEnd;
    [~,IA,~] = intersect(CSPtsRtC0,PtsOnTT,'rows','stable');
    WeightPtsTT = Weights(IA,:);
    WeightPtsTT = WeightPtsTT/sum(WeightPtsTT);
    PtMeanTT = sum([WeightPtsTT WeightPtsTT WeightPtsTT].*PtsOnTT,1);
    
    IDX = knnsearch(CSPtsRtC0,PtMeanTT);
    PtTT = CSPtsRtC0(IDX,:);
    UTT = PtTT'/norm(PtTT);
    
    % Identify the medial third of the Tibial Tuberosity :
    % Make it lateral to medial like for the Tibia CS
    VTT = [UTT(2) -UTT(1) 0]'; VTT = sign(VTT(2))*VTT;
    [~,IMedialPt] = max(PtsOnTT*VTT);
    [~,ILateralPt] = min(PtsOnTT*VTT);
    
    MedialThirdOfTT = 1/3*PtsOnTT(ILateralPt,:) + 2/3*PtsOnTT(IMedialPt,:);
    
    LineEndingOfTT = PtsOnTT(IMedialPt,:)-PtsOnTT(ILateralPt,:);
    VLineEndingOfTT = LineEndingOfTT'/norm(LineEndingOfTT);
    
    [~,IDMedThird] = min(abs(bsxfun(@minus,PtsOnTT,MedialThirdOfTT)*VLineEndingOfTT));
    
    PtMedialThirdOfTT = PtsOnTT(IDMedThird,:);
    
    
    if plots ==1
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
    PtsTT(j,:) = transpose(CS.V*(PtTT + Centroid)') ;
    
end

IDPtMedialThird = knnsearch(ProxTib.Points,mean(PtsMedThird));
PtMedialThirdOfTT = ProxTib.Points(IDPtMedialThird,:);

if plots ==1
    figure()
    trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',1,'edgecolor','none'); % 0.8,0.8,0.85
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