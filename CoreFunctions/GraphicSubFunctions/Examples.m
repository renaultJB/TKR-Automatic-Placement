
%% Example of Plot

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
    
    
    
    figure()
    trisurf(ProxTib,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.5,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    axis equal
    light('Position',CS.Origin' + 300*CS.Y + 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 200*CS.Y - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 50*CS.Y + 50*CS.X - 500*CS.Z,'Style','local')
    hold on
    grid off
    lighting gouraud
    pl3t(BoundaryStemTip,'b-','linewidth',3)
    pl3t(Boundary_xp,'r-','linewidth',3)
    plotDot( StemTip_CT', 'g', 1.5 )
    plotDot( CDiaphysisStemTip_CT, 'r', 1.5 )
    plotDot(PtMedialThirdOfTT)
    plotDot( ProthOrig, 'm', 1.5 )
    trisurf(ProsthesisEnd,'Facecolor','g','FaceAlpha',1,'edgecolor','none');
    
    
        
    
    figure()
    trisurf(Prosthesis,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.5,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    axis equal
    light('Position',[50 50 200],'Style','local')
    light('Position',[-50 -50 -200],'Style','local')
    light('Position',[-50 100 -200],'Style','local')
    hold on
    grid off
    lighting gouraud
    pl3t(Boundary_xp_inRxp,'b-','linewidth',3)
    pl3t(Boundary_ProsthesisTP,'r-')
    plotDot(TTproj)