function PlotPosOptim2D( x, Prosthesis, Boundary_xp, ProxTib, GS, GS_TTA,...
    PtMiddleOfTT, Oxp, PtMedialThirdOfTT, R_xp, CS, history )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% mean(ProxTibCuttedCS.Points)

if nargin < 12
    history = x;
end

if size(history,1) > 5
    v = VideoWriter('newfile.mp4','MPEG-4');
    v.FrameRate = 8;
    open(v)
end

%% Prepare some data

Nxp = R_xp(:,3);

%

Curve_ttam = TriPlanIntersect(ProxTib,CS.Z,PtMiddleOfTT);

Curve_mtttam = TriPlanIntersect(ProxTib,CS.Z,-CS.AltStartEpiph);


Curve_TP = TriPlanIntersect(ProxTib,CS.nTP,CS.Origin-5*CS.nTP');


%%

DeltaTP = 1/4 * (CS.Origin*CS.Z-CS.AltStartEpiph);

ProxTib = triangulation(ProxTib.ConnectivityList,...
    transpose(CS.V'*bsxfun(@minus,ProxTib.Points,CS.Origin)'));

TP_Pts = ProxTib.Points(ProxTib.Points(:,3)>-DeltaTP,:); %-DeltaTP

k = convhull(TP_Pts(:,1),TP_Pts(:,2));

TP_CH_Pts = [TP_Pts(k,1),TP_Pts(k,2),zeros(length(k),1)];



% Get implant tray contour
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
% Simplify Implant contour
Boundary_ProsthesisTP = [CurvesProsthesisTP.Pts(1:6:end-1,:) ; CurvesProsthesisTP.Pts(end,:)];

%%
figure(1)
figure('units','normalized','outerposition',[0 0 1 1])
view([0,90])
grid off

for j = 1 : size(history,1)
    x = history(j,:);
    
    %% Plot Fixed bouindaries
    Boundary_xp_inRt = CS.V'*bsxfun(@minus,Boundary_xp,CS.Origin)';
    hold off
    pl3t(Boundary_xp_inRt,'g-','linewidth',4)
    hold on;
    
%     pl3t(TP_CH_Pts,'k*-')
    
    % Plot Implant in place
    gamma = deg2rad(x(3));
    R = [cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0; 0 0 1];
    
    %% Tibial implant Contour
    ProsthContourTR_tmp = R*Boundary_ProsthesisTP';
    ProsthContourTR_optim = bsxfun(@plus,ProsthContourTR_tmp',[x(1) x(2) 0]);
    ProsthContourTR_optim = R_xp*ProsthContourTR_optim';
    ProsthContourTR_optim = bsxfun(@plus,ProsthContourTR_optim',Oxp);
    IT_Contour_optim_inRt = CS.V'*bsxfun(@minus,ProsthContourTR_optim,CS.Origin)';
    
    pl3t(IT_Contour_optim_inRt,'r-','linewidth',2)
    
    O_it_inRt = CS.V'*(Oxp'+R_xp*[x(1) x(2) 0]'-CS.Origin');
    % Above Tibial tray
    for dZ = [ 0.15, 0.40, 0.80, 1.50, 2.5]
        ProsthesisShape2 = TriPlanIntersect(Prosthesis,[10^-6; 10^-6; 1],-dZ);
        
        for i=1:length(ProsthesisShape2)
            ProsthesisShape2(i).Pts = R*ProsthesisShape2(i).Pts';
            ProsthesisShape2(i).Pts = bsxfun(@plus,ProsthesisShape2(i).Pts',[x(1) x(2) 0]);
            
            ProsthesisShape2(i).Pts = R_xp*ProsthesisShape2(i).Pts';
            ProsthesisShape2(i).Pts = bsxfun(@plus,ProsthesisShape2(i).Pts',Oxp);
            ProsthesisShape2(i).Pts = CS.V'*bsxfun(@minus,ProsthesisShape2(i).Pts,CS.Origin)';
            
            pl3t(ProsthesisShape2(i).Pts,'r-')
        end
    end
    
    %% Cut Orientaton
    Oxp_inRt = CS.V'*(Oxp-CS.Origin)';
    
    Nxp_inRt = CS.V'*Nxp;
    plotArrow(Nxp_inRt,0.25,Oxp_inRt,40,1,'r')
    plotArrow([0;0;1],0.25,Oxp_inRt,40,1,'k')
    
    %% Plot TP contour account for tibia deformation
    Curve_TP_inRt = CS.V'*bsxfun(@minus,Curve_TP(1).Pts,CS.Origin)';
%     pl3t(Curve_TP_inRt,'g--','linewidth',1.5)
    
    IDPtMostAntofPost = knnsearch(Curve_TP_inRt',[0,CS.ML_Width_xp/12,0]);
    Pt_Akagi_PostTib = Curve_TP_inRt(:,IDPtMostAntofPost)';
    
    dRight = Pt_Akagi_PostTib(1)-min(Curve_TP_inRt(1,:));
    dLeft = max(Curve_TP_inRt(1,:))-Pt_Akagi_PostTib(1);
    d = min(dRight,dLeft);
    Pt_Akagi_PostTib = Pt_Akagi_PostTib + [sign(Pt_Akagi_PostTib(1))*d/2,0,0];
    
    PT_Insall_PostTib = Pt_Akagi_PostTib;
    plotDot(Pt_Akagi_PostTib,'c',1)
    
    %% Plot countour at TTA
    GS_inRt = CS.V'*(GS'-CS.Origin)';
    
    %%------------------------
%     GS_inRt = [Oxp_inRt(1:2);0]
    %%---------------------------
    
    
    GS_TTA_inRt = CS.V'*(GS_TTA'-CS.Origin)';
    PtMiddleOfTT_inRt = CS.V'*(PtMiddleOfTT-CS.Origin)';
    PtMedialThirdOfTT_inRt = CS.V'*(PtMedialThirdOfTT-CS.Origin)';
    
    
    
    Curve_ttamPts = Curve_ttam(1).Pts;
    Curve_ttam_inRt = CS.V'*bsxfun(@minus,Curve_ttamPts,CS.Origin)';
    
    Curve_mtttamPts = Curve_mtttam(1).Pts;
    Curve_mtttam_inRt = CS.V'*bsxfun(@minus,Curve_mtttamPts,CS.Origin)';
    % pl3t(Curve_ttam_inRt,'b-')
    plotDot(PtMiddleOfTT_inRt','b',1.5)
    
%     plotDot(PtMedialThirdOfTT_inRt','g',1)
    pl3t(GS_inRt,'ko')
    pl3t(GS_TTA_inRt,'ko')
    % Keep only exterior part of the section at the TTA level
    
    d = p_poly_dist(Curve_ttam_inRt(1,:), Curve_ttam_inRt(2,:),...
        Boundary_xp_inRt(1,:), Boundary_xp_inRt(2,:));
    pl3t(Curve_ttam_inRt,'g--')
    pl3t(Curve_ttam_inRt(:,d>0),'g.')
    
    d = p_poly_dist(Curve_mtttam_inRt(1,:), Curve_mtttam_inRt(2,:),...
        Boundary_xp_inRt(1,:), Boundary_xp_inRt(2,:));
    
    Curve_mtttam_inRt_out = Curve_mtttam_inRt(:,d>0);
    [~,Id_Akagi] = max( Curve_mtttam_inRt_out(2,:));
    [~,Id_InsallLAt] = min( Curve_mtttam_inRt_out(2,:));
    Pt_TT_Akagi = Curve_mtttam_inRt_out(:,Id_Akagi)';
    plotDot(Pt_TT_Akagi,'c',1)
    
    Pt_TT_Insall_Med = Pt_TT_Akagi;
    Pt_TT_Insall_Lat = Curve_mtttam_inRt_out(:,Id_InsallLAt)';
    Pt_TT_Insall = 2/3*Pt_TT_Insall_Med + 1/3*Pt_TT_Insall_Lat;
    plotDot(Pt_TT_Insall,'m',1)

    pl3t(Curve_mtttam_inRt,'g--')
%     pl3t(Curve_mtttam_inRt(:,d>0),'g.')
    
    %% Plot Lines of GS to TTA
    
    U_TTA_inRt = normalizeV(PtMiddleOfTT_inRt-GS_TTA_inRt);
    U_MTTTA_inRt = normalizeV(PtMedialThirdOfTT_inRt-GS_TTA_inRt);
    

    LineTT_Akagi = [Pt_TT_Akagi',Pt_Akagi_PostTib'];
    LineTT_Insall = [Pt_TT_Insall',PT_Insall_PostTib'];
    pl3t(LineTT_Akagi,'c--','linewidth',1.5)
    pl3t(LineTT_Insall,'m--','linewidth',1.5)
    
    LineTTA_inRT = [GS_TTA_inRt,PtMiddleOfTT_inRt];
    LineMTTTA_inRT = [GS_TTA_inRt,PtMedialThirdOfTT_inRt];
    
    pl3t(LineTTA_inRT,'b-','linewidth',2)
%     pl3t(LineMTTTA_inRT,'k-','linewidth',2)
    
    X_it = CS.V'*R_xp*R(:,1);
    X_it = sign(X_it'*U_TTA_inRt)*X_it;
    
    X_it_proj = normalizeV([X_it(1:2);0]);
    
    Line_it = [O_it_inRt,O_it_inRt+45*X_it_proj];
    Line_it_GS = [GS_inRt,GS_inRt+40*X_it_proj];
    pl3t(Line_it,'r-','linewidth',1)
    pl3t(Line_it_GS,'r-','linewidth',2)
    
    

    view([0,90])
    axis off
    if size(history,1) > 5
        writeVideo(v,getframe(gcf));
    end
    
end

if size(history,1) > 5
    for i=1:round(length(history)/5)
        writeVideo(v,getframe(gcf));
    end
    
    close(v)
    
end


end

