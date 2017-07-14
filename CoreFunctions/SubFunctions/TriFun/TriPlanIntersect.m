function [ Curves , Area , InterfaceTri ] = TriPlanIntersect( Tr, n , d )
%TriPlanIntersect Intersection Triangle Plane in 3D
%   Detailed explanation goes here


Nodes = Tr.Points;
Elmts = Tr.ConnectivityList;

% Normalize vector and compute d altitude for (0,0,z) if a point is given
% for input argument d



n=n/norm(n);
d=double(d);
% If d is a point on the plan and not the altitude at (0,0,z)
if length(d)>2
    if size(d,2)==1
        d = d'*n;
    elseif size(d,1)==1
        d = d*n;
    else
        sprintf('ERROR')
    end
    
end

% Particularity to deal with round off errors
antiRoundOff = 0;
%unexplained bug
if abs(d)<1
    antiRoundOff = 1 ;
end

d = d + 10*antiRoundOff ;
Nodes = bsxfun(@plus, Nodes , 10*antiRoundOff*n') ;

NodesDist = Nodes*n-d;
Elmts_Status = sum(sign(NodesDist(Elmts)),2);

IndexInterfaceTri = abs(Elmts_Status) < 3;

if ~logical(mean(IndexInterfaceTri))
    Area = 0;
    InterfaceTri = [];
    Curves = struct();
    Curves(1).NodesID = [];
    Curves(1).Pts = [];
    Curves(1).Text = 'No Intersection';
else
    InterfaceTri = Elmts(IndexInterfaceTri,:);
    % Create intersection Points
    PtsInter = zeros(length(InterfaceTri)*2,3);
    Segments = zeros(length(InterfaceTri),2);
    j=1;
    kp=0;
    for k = 1: length(InterfaceTri)
        da = NodesDist(InterfaceTri(k,1));
        db = NodesDist(InterfaceTri(k,2));
        dc = NodesDist(InterfaceTri(k,3));
        A = Nodes(InterfaceTri(k,1),:);
        B = Nodes(InterfaceTri(k,2),:);
        C = Nodes(InterfaceTri(k,3),:);
        % on allow to deal with elements with only  a node on the plan or
        % all node on plan
        on = 0;
        
        if da*db < 0
            PtsInter(j,:) = A + da/(da-db)*(B-A);
            j=j+1;
            on = on + 1;
        elseif da == 0
            PtsInter(j,:) = A;
            j=j+1;
            on = on + 1;
        end
        if db*dc < 0
            PtsInter(j,:) = B + db/(db-dc)*(C-B);
            j=j+1;
            on = on + 1;
        elseif db == 0
            PtsInter(j,:) = B;
            j=j+1;
            on = on + 1;
        end
        if da*dc < 0
            PtsInter(j,:) = C + dc/(dc-da)*(A-C);
            j=j+1;
            on = on + 1;
        elseif dc == 0
            PtsInter(j,:) = C;
            j=j+1;
            on = on + 1;
        end
        
        if on == 1
            PtsInter(j,:) = [] ;
            j=j-1;
            PtsInter(end,:) = [];
            Segments(end,:) = [];
            kp=kp+1;
        elseif on == 3
            PtsInter(j-2:j,:) = [] ;
            j=j-3;
            kp=kp+3;
        else
            Segments(k-kp,:) = [ j-2 j-1] ;
        end
    end
    
    [PtsOut, ~ , indexn] =  uniquetol(PtsInter , eps('double') , 'ByRows',true);
    
    Segments = indexn(Segments);
    
    
    if numel(PtsOut) ~= numel(PtsInter)/2
        Segments(Segments(:,1) == Segments(:,2),:) = []  ;
    end
    
    j=1;
    Curves=struct();
    i=1;
    while ~isempty(Segments)
        Curves(i).NodesID = zeros(length(PtsInter),1);
        Curves(i).NodesID(j)=Segments(1,1);
        Curves(i).NodesID(j+1)=Segments(1,2);
        Segments(1,:)=[];
        j=j+1;
        [Is,Js] = ind2sub(size(Segments),find(Segments(:) == Curves(i).NodesID(j)));
        Nk = Segments(Is,round(Js+2*(1.5-Js)));
        Segments(Is,:)=[];
        j=j+1;
        while ~isempty(Nk)
            
            Curves(i).NodesID(j) = Nk;
            [Is,Js] = ind2sub(size(Segments),find(Segments(:) == Curves(i).NodesID(j)));
            Nk = Segments(Is,round(Js+2*(1.5-Js)));
            Segments(Is,:)=[];
            j=j+1;
        end
        i=i+1;
    end
    
    Area = 0;
    MatrixArea = zeros(length(Curves));
    for i = 1 : length(Curves)
        Curves(i).NodesID(Curves(i).NodesID==0) = [];
        Curves(i).Pts = bsxfun(@minus, PtsOut(Curves(i).NodesID,:) , 10*antiRoundOff*n') ;
        % Replace the close curve in coordinate system where X, Y or Z is 0, in
        % % order to use polyarea
        [V,~]= eig(cov(Curves(i).Pts));
        CloseCurveinRplanar = V'*Curves(i).Pts';
        % First principal inertial vector is normal to plan
        Area = Area + polyarea(CloseCurveinRplanar(2,:),CloseCurveinRplanar(3,:));
        MatrixArea(i,i) = polyarea(CloseCurveinRplanar(2,:),CloseCurveinRplanar(3,:));
    end
    
    for i = 1 : length(Curves)
        for j =  1 : length(Curves)
            if i~=j
                Curves(i).NodesID(Curves(i).NodesID==0) = [];
                Curves(i).Pts = bsxfun(@minus, PtsOut(Curves(i).NodesID,:) , 10*antiRoundOff*n') ;
                % Replace the close curve in coordinate system where X, Y or Z is 0, in
                % % order to use polyarea
                [V,~]= eig(cov(Curves(i).Pts));
                CloseCurveinRplanar1 = V'*Curves(i).Pts';
                
                Curves(j).NodesID(Curves(j).NodesID==0) = [];
                Curves(j).Pts = bsxfun(@minus, PtsOut(Curves(j).NodesID,:) , 10*antiRoundOff*n') ;
                % Replace the close curve in coordinate system where X, Y or Z is 0, in
                % order to use polyarea
                CloseCurveinRplanar2 = V'*Curves(j).Pts';
                
                if inpolygon(CloseCurveinRplanar2(2,1),CloseCurveinRplanar2(3,1),...
                        CloseCurveinRplanar1(2,:),CloseCurveinRplanar1(3,:))>0
                    MatrixArea(i,j) = -2*polyarea(CloseCurveinRplanar2(2,:),CloseCurveinRplanar2(3,:));
                end
            end
        end
    end
    
    
    Area = sum(sum(MatrixArea));
    
    
end




%% Inset of the extorior polygon

end

