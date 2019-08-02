function [ x,fval,history ] = optimC_PlacementTI_xp( x0,Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA )
%Optimisation of the tibial implant placement in the resection face
%  Try to follow the constrained presented by Berger et al. 1998


%     f = @(x)CenteringStemTipGoal(x , Boundary_xp_inRxp , Boundary_ProsthesisTP , StemTip);  % C or CDiaphysisStemTip
%     fcon = @(x)CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP , TTproj);
%     x = fmincon(f,x0,A,b,Aeq,beq,lb,ub,fcon);

history = [];

f = @(x)CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP , CS, R_xp, theta_TTA);
options = optimoptions(@fminunc,'Algorithm','quasi-newton','MaxFunctionEvaluations',500,'OutputFcn', @myoutput);
[x,fval] = fminunc(f,x0,options);

    function stop = myoutput(x,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
            history = [history; x];
        end
    end

    function [ C,Ceq ] = CoverageCost(x, Boundary_xp_inRxp, Boundary_ProsthesisTP, CS, R_xp, theta_TTA)
        %CoverageCost try to get a cost function for respecting the no
        %overhang policy and penalise bad alignement
        %   
        
        x(3) = x(3)*pi/180 ; % Just to have x(1),x(2),x(3) with the same order of magnitude
        Boundary_xp_inRxp(:,3) = [];
        Boundary_ProsthesisTP(:,3) = [];
        
        ProsthOrig = [x(1) x(2)];
        
        %% Displace implant
        % 1st rotation around original axis
        R = [cos(x(3)) -sin(x(3));sin(x(3)) cos(x(3))];
        ProsthContourTR_tmp = R*Boundary_ProsthesisTP';
        % 2nd translate origin
        ProsthContourTR = bsxfun(@plus,ProsthContourTR_tmp',ProsthOrig);
        
        
        
        
        
        %% Compute placement criteria
        % Overhang   
        d = p_poly_dist(ProsthContourTR(:,1), ProsthContourTR(:,2), Boundary_xp_inRxp(:,1), Boundary_xp_inRxp(:,2));
        
        % Rotational Alignement
        U_it = [cos(x(3)) ; sin(x(3)) ; 0];
        U_t = normalizeV(CS.Paxial*(R_xp*U_it));
        theta_it = rad2deg(acos(CS.Y'*U_t));
        deltaTheta = (theta_TTA - 10) - theta_it; % Rotational error in degree 18° is the physiological value
        
        %% Compute cost function
        % Criterium 1 : No overhang of the implant
        C1 = 30 * mean(exp(d).^2);  % Cost function , d+0.25 for penalty if the prosthesis is too close of the edges
        
        % Criterium 2 : Implant oriented towards the Tibial tuberosity
        C2 = (4*abs(deltaTheta) + (deltaTheta)^2)/5;
        
        % Overall, sum of the two criteria
        C = C1 + C2;
        
        Ceq = -1;
    end

end
