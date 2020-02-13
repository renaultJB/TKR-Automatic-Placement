function [ptOnStemTip] = getPtOnStemTip(implant,nXP)
%Get a point on the plane at the stem tip of the implant
%   XP : Cut Plan

minDistToXP = min(implant.incenter*nXP);

IdElmtsOk = find(implant.faceNormal*nXP < -0.985 &...
    implant.incenter*nXP < (minDistToXP + 5));



StemTip = TriReduceMesh(implant,IdElmtsOk);

% figure()
% trisurf(StemTip);
% axis equal

% 2nd iteration
[x0,nST] = lsplane(StemTip.incenter);
nST = sign(nST'*nXP)*nST;
acosd(nST'*nXP);

IdElmtsOk = find(StemTip.faceNormal*nST < -0.995 &...
    StemTip.incenter*nST <= mean(StemTip.incenter)*nST);
StemTip = TriReduceMesh(StemTip,IdElmtsOk);

IDX = knnsearch(StemTip.incenter,mean(StemTip.incenter));

Pts = StemTip.incenter;



ptOnStemTip = Pts(IDX,:);

% figure()
% trisurf(StemTip);
% plotDot(ptOnStemTip,'k',2)
% axis equal


end

