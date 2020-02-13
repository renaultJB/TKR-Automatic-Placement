function [ x,fval,exitflag,output,history ] = optimC_PlacementTI_xp( x0, A,b,Aeq,beq,lb,ub, Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA )
%Optimisation of the tibial implant placement in the resection face
%  Try to follow the constrained presented by Berger et al. 1998


%     f = @(x)CenteringStemTipGoal(x , Boundary_xp_inRxp , Boundary_ProsthesisTP , StemTip);  % C or CDiaphysisStemTip
%     fcon = @(x)CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP , TTproj);
%     x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,fcon);

history = [];

f = @(x)AxialRotation(x , CS, R_xp, theta_TTA);
fcon = @(x)CoverageConstraint(x, Boundary_xp_inRxp, Boundary_ProsthesisTP);
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxFunctionEvaluations',5000,'OutputFcn', @myoutput);
[x,fval,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,fcon,options);

    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
            history = [history; x];
        end
    end

    function [ C,Ceq ] = AxialRotation(x, CS, R_xp, theta_TTA)
        %CoverageCost try to get a cost function for respecting the no
        %overhang policy and penalise bad alignement
        %   
        
        x(3) = x(3)*pi/180 ; % Just to have x(1),x(2),x(3) with the same order of magnitude

        %% Compute orientation criteria
        % Rotational Alignement
        U_it = [cos(x(3)) ; sin(x(3)) ; 0];
        U_t = normalizeV(CS.Paxial*(R_xp*U_it));
        theta_it = rad2deg(acos(CS.Y'*U_t));
        deltaTheta = theta_TTA - theta_it; % Rotational error in degree 18° is the physiological value
        
        %% Compute objective function
        % Criterium 2 : Implant oriented towards the Tibial tuberosity
        C = (4*abs(deltaTheta) + (deltaTheta)^2)/5;
        
        Ceq = 0;
    end




    function [ C,Ceq ] = CoverageConstraint(x, Boundary_xp_inRxp, Boundary_ProsthesisTP)
        
        x(3) = x(3)*pi/180 ; % Just to have x(1),x(2),x(3) with the same order of magnitude
        Boundary_xp_inRxp(:,3) = [];
        Boundary_ProsthesisTP(:,3) = [];
        ProsthOrig = [x(1) x(2)];
        
        % Smooth high curvature in the cut plan R around 5 to 8 mm
        SmoothingRadius = 0.1*range(Boundary_xp_inRxp(:,2));
        Boundary_xp_inRxp2 = [Boundary_xp_inRxp;...
            0.8*Boundary_xp_inRxp;...
            0.6*Boundary_xp_inRxp;...
            0.4*Boundary_xp_inRxp;...
            0.2*Boundary_xp_inRxp];

        shp = alphaShape(Boundary_xp_inRxp2(:,1), Boundary_xp_inRxp2(:,2),...
                         SmoothingRadius,'HoleThreshold',10^5);
        BF = boundaryFacets(shp);
        ID=BF(:,1);
        
        %% Displace implant
        % 1st rotation around original axis
        R = [cos(x(3)) -sin(x(3));sin(x(3)) cos(x(3))];
        
        ProsthContourTR_tmp = R*Boundary_ProsthesisTP';
        % 2nd translate origin
        ProsthContourTR = bsxfun(@plus,ProsthContourTR_tmp',ProsthOrig);
        
        %% Compute placement criteria
        % Overhang
        d = p_poly_dist(ProsthContourTR(:,1), ProsthContourTR(:,2), Boundary_xp_inRxp2(ID,1), Boundary_xp_inRxp2(ID,2));
        
        %         figure(100)
        %         plot( Boundary_xp_inRxp2(ID,1), Boundary_xp_inRxp2(ID,2),'k*-')
        %         axis equal
        %         view([0,90])
        
        
        %% Compute cost function
        % Criterium : No overhang of the implant
        C = max(exp(d+0.25)-1);  % Cost function , d+0.25 for penalty if the prosthesis is too close of the edges
        %         C = max(d+0.25);
%         % Align medioLateral, good repartition of implant distance to the
%         % medial and lateral border of the tibia
%         delta_Med = max(Boundary_xp_inRxp(:,2))-max(ProsthContourTR(:,2));
%         delta_Lat = min(ProsthContourTR(:,2))-min(Boundary_xp_inRxp(:,2));
%         RatioML = abs(delta_Med/(delta_Med+delta_Lat)-1/2)
%         
%         % Could also use " 'ConstraintTolerance',1/6 " in option
%         if RatioML<1/6
%             Ceq = 0;
%         else
%             Ceq = RatioML;
%         end
%         
        Ceq = [];
        
        
    end

end
