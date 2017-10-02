function [ CentralityScore, T ] = CentralityScore(ProxTib, Prosth0, Prosth, StemTip, LegSide)
% CENTRALITYSCORE evaluate the centrality of the stem tip relative to the 
% Tibial Bone

%% 1st Identify the ring of the the stem tip

IdElmtsOk = find(Prosth0.faceNormal*[0 0 -LegSide]'< 0.985 &...
    Prosth0.faceNormal*[0 0 -LegSide]' > 0.18 &...
    Prosth0.incenter*[0 0 1]' < 0.951 * StemTip(3));

Ring = TriReduceMesh(Prosth,IdElmtsOk);



%% 2nd identify a layer of the proxTibia at the same altitude
RingPts = Ring.Points;

[n,d] = LS_Plan(RingPts);

IdElmtsOk =  find(abs(ProxTib.incenter*n+(d+1)) <3);

ProxTibLayer = TriReduceMesh(ProxTib,IdElmtsOk);

ProxTibLayerPts = ProxTibLayer.Points;






%% Compute the paired points closest distance



%% 
[IDX,D] = knnsearch(RingPts,ProxTibLayerPts);

RingPtsOk = RingPts(IDX,:);


CentralityScore.CV = 1-std(D)/mean(D);
CentralityScore.MinMax = min(D)/max(D);
CentralityScore.MinMean = min(D)/mean(D);

%% Create Table

RowNames = {'Xring';'Yring';'Zring';'Xtib';'Ytib';'Ztib';'Distance'};
[Xring,Yring,Zring,Xtib,Ytib,Ztib,Distance] = deal(RingPtsOk(:,1),RingPtsOk(:,2),RingPtsOk(:,3),...
    ProxTibLayerPts(:,1),ProxTibLayerPts(:,2),ProxTibLayerPts(:,3),D);

T = table(Xring,Yring,Zring,Xtib,Ytib,Ztib,Distance,'VariableNames',RowNames);


 
% Plot
% figure()
% trisurf(ProxTibLayer,'Facecolor',[0.65    0.65    0.6290],'FaceAlpha',0.5,'edgecolor','none'); % 0.8,0.8,0.85
% hold on
% axis equal
% trisurf(Ring,'Facecolor','r','FaceAlpha',0.5,'edgecolor','none'); % 0.8,0.8,0.85
% light('Position',mean(Ring.Points) + [0 0 300],'Style','local')
% light('Position',mean(Ring.Points) + [100 -100 100],'Style','local')
% light('Position',mean(Ring.Points) + [-100 100 -100],'Style','local')
% light('Position',mean(Ring.Points) + [-100 -100 -100],'Style','local')
% hold on
% grid off
% lighting gouraud
% 
% hold on 
% for i = 1:length(Xring)
%     plot3([Xring(i,:) Xtib(i,:)],[Yring(i,:) Ytib(i,:)],[Zring(i,:) Ztib(i,:)],'b-')
% end


end