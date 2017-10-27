function PlotPosOptim( ProxTib, Prosthesis0, history, Start_Point, Oxp, U_xp, V_xp, Nxp, R_xp, LegSide, d_xp, CS, PtMedialThirdOfTT)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Remove ProxTib elements that are above the cut line
ElemtsOk = find(ProxTib.incenter*Nxp < Oxp*Nxp);
ProxTibCutted = TriReduceMesh(ProxTib,ElemtsOk);
[ProxTibCutted, HolePtsProj]  = TriFillPlanarHoles( ProxTibCutted, 1);

v = VideoWriter('newfile.mp4','MPEG-4');
v.FrameRate = 8;
open(v)

figure(1)

figure('units','normalized','outerposition',[0 0 1 1])
% trisurf(ProxTibCutted,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',00,'edgecolor','none'); % 0.8,0.8,0.85
% hold on
% axis equal
% light('Position',CS.Origin' + 300*CS.Y + 200*CS.X,'Style','local')
% light('Position',CS.Origin' + 200*CS.Y - 200*CS.X,'Style','local')
% light('Position',CS.Origin' + 50*CS.Y + 50*CS.X - 500*CS.Z,'Style','local')
% plotDot( PtMedialThirdOfTT, 'r', 2 )

for i = 1 : length(history)
    x = history(i,:);
    ProthOrig = Start_Point + x(1)*U_xp' + x(2)*V_xp';
    Rp = rot(Nxp,x(3));
    
    PtsProsth0 = Prosthesis0.Points;
    PtsProsth0(:,4) = ones(length(PtsProsth0),1);

    T = zeros(4,4); T(1:3,1:3) = Rp*R_xp*[0 LegSide 0 ; LegSide 0 0; 0 0 -1]; %[0 LegSide 0 ; 1 0 0; 0 0 -1]
    T(:,4)=[ProthOrig';1];
    
    
    %% Move Prosthesis points
    PtsProsthEnd = transpose(T*PtsProsth0');
    PtsProsthEnd(:,4)=[];
    
    ProsthesisEnd = triangulation(Prosthesis0.ConnectivityList,PtsProsthEnd);
    
    
    %% Plot
    subplot(1,5,1:3)
    hold off
    trisurf(ProxTibCutted,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.75,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    pl3t(HolePtsProj,'b.')
    axis equal
    light('Position',CS.Origin' + 1000*CS.Z + 200*CS.X,'Style','local')
%     light('Position',CS.Origin' + 200*CS.Y - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 550*CS.Y - 500*CS.X + 800*CS.Z,'Style','local')
    plotDot( PtMedialThirdOfTT, 'r', 2 )
    trisurf(ProsthesisEnd,'Facecolor','g','FaceAlpha',1,'edgecolor','none');
    plotDot( PtMedialThirdOfTT, 'r', 2 )
    hold on
    grid off
    axis off
%     lighting gouraud
    view([180 90])
    
    subplot(1,5,4)
        hold off
    trisurf(ProxTibCutted,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.6,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    pl3t(HolePtsProj,'b.')
    axis equal
    light('Position',CS.Origin' + 1000*CS.Z + 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 200*CS.Y - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 550*CS.Y - 500*CS.X + 800*CS.Z,'Style','local')
    plotDot( PtMedialThirdOfTT, 'r', 2 )
    trisurf(ProsthesisEnd,'Facecolor','g','FaceAlpha',1,'edgecolor','none');
    plotDot( PtMedialThirdOfTT, 'r', 2 )
    hold on
    grid off
    axis off
    lighting gouraud
    view([0 0])
    
        subplot(1,5,5)
        hold off
    trisurf(ProxTibCutted,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.6,'edgecolor','none'); % 0.8,0.8,0.85
    hold on
    pl3t(HolePtsProj,'b.')
    axis equal
    light('Position',CS.Origin' + 1000*CS.Z + 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 200*CS.Y - 200*CS.X,'Style','local')
    light('Position',CS.Origin' + 550*CS.Y - 500*CS.X + 800*CS.Z,'Style','local')
    plotDot( PtMedialThirdOfTT, 'r', 2 )
    trisurf(ProsthesisEnd,'Facecolor','g','FaceAlpha',1,'edgecolor','none');
    plotDot( PtMedialThirdOfTT, 'r', 2 )
    hold on
    grid off
    axis off
    lighting gouraud
    view([90 0])
    
    writeVideo(v,getframe(gcf));    
    
end

close(v)

end

