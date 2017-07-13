function [ TRout ] = TriDifferenceMesh( TR1 , TR2 )
%Boolean difference between original mesh TR1 and another mesh TR2
% /!\ delete all elements in TR1 that contains a node in TR2

[~, ia , ~] = intersect(TR1.Points,TR2.Points,'rows','stable');

Elmts2Delete = TR1.vertexAttachments(ia)';
ElmtsAll = ones(length(TR1.ConnectivityList),1);

Elmts2Keep = ElmtsAll;
Elmts2Keep(unique(horzcat(Elmts2Delete{:}))) = 0;
Elmts2KeepID = find(Elmts2Keep);
TRout = TriReduceMesh( TR1, Elmts2KeepID);


end

