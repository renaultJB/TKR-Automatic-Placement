%%-------------------------------------------------------------------------  
% GOAL :
% Read the results and compute the minimal bounding sphere radii
% of the origins of each algorithm, for the patella,femur and tibia.

% Output ? 3 tables :
% - mBSF.txt
% - mBST.txt
% - mBSP.txt
% Tables will be subsequently exploited in R
%% ------------------------------------------------------------------------
addpath(strcat(pwd,'\MinBoundSphere'));
clearvars

% Load results generated by the ACS_LL.m script
load('Results.mat')
N = length(Results);

% Initialize the Bounding Spheres radii, N by P matrices with, N?#subjects, P?#algorithms 
BSF = zeros(N,5); % F for Femur 5 algo
BST = zeros(N,5); % T for Tibia 5 algo
BSP = zeros(N,4); % P for Patella 4 algo

for i = 1 : N
    %% Femur
    
    Ccyl =  [Results(i, 1).Fem.CenterKnee; Results(i, 2).Fem.CenterKnee; Results(i, 3).Fem.CenterKnee];
    Csph =  [Results(i, 1).Fem.CenterKneeSph'; Results(i, 2).Fem.CenterKneeSph'; Results(i, 3).Fem.CenterKneeSph'];
    Celpsd =  [Results(i, 1).Fem.CenterKneeElpsd'; Results(i, 2).Fem.CenterKneeElpsd'; Results(i, 3).Fem.CenterKneeElpsd'];
    CKai =  [Results(i, 1).Fem.CenterKneeKai; Results(i, 2).Fem.CenterKneeKai; Results(i, 3).Fem.CenterKneeKai];
    CMir = [Results(i, 1).Fem.CenterKnee_Miranda; Results(i, 2).Fem.CenterKnee_Miranda; Results(i, 3).Fem.CenterKnee_Miranda];
    
    % Because the ApproxMinBoundSphereND function need 4 non-coplanar points
    % a 4th point is added, it is the barycenter of the 3 origin points
    % randomly perturbated in space :
    Ccyl = [Ccyl; mean(Ccyl)+0.005*rand(1,3)];
    Csph = [Csph; mean(Csph)+0.005*rand(1,3)];
    Celpsd = [Celpsd; mean(Celpsd)+0.005*rand(1,3)];
    CKai = [CKai; mean(CKai)+0.005*rand(1,3)];
    CMir = [CMir; mean(CMir)+0.005*rand(1,3)];
    
    
    % Get the minimal bounding sphere radii of each algorithm:
    [BSF(i,1),~] = ApproxMinBoundSphereND(Ccyl);
    [BSF(i,2),~] = ApproxMinBoundSphereND(Csph);
    [BSF(i,3),~] = ApproxMinBoundSphereND(Celpsd);
    [BSF(i,4),~] = ApproxMinBoundSphereND(CKai);
    [BSF(i,5),~] = ApproxMinBoundSphereND(CMir);
    
    clearvars CKai CMir
    
    %% Tibia
    
    C1 =  [Results(i, 1).Tib.tech1.CenterKnee; Results(i, 2).Tib.tech1.CenterKnee; Results(i, 3).Tib.tech1.CenterKnee];
    C2 =  [Results(i, 1).Tib.tech2.CenterKnee; Results(i, 2).Tib.tech2.CenterKnee; Results(i, 3).Tib.tech2.CenterKnee];
    C3 =  [Results(i, 1).Tib.tech3.CenterKnee; Results(i, 2).Tib.tech3.CenterKnee; Results(i, 3).Tib.tech3.CenterKnee];
    CKai =  [Results(i, 1).Tib.Kai.CenterKnee; Results(i, 2).Tib.Kai.CenterKnee; Results(i, 3).Tib.Kai.CenterKnee];
    CMir =  [Results(i, 1).Tib.Miranda.CenterKnee; Results(i, 2).Tib.Miranda.CenterKnee; Results(i, 3).Tib.Miranda.CenterKnee];
    
    % Because the ApproxMinBoundSphereND function need 4 non-coplanar points
    % a 4th point is added, it is the barycenter of the 3 origin points
    % randomly perturbated in space :
    C1 = [C1;mean(C1)+0.005*rand(1,3)];
    C2 = [C2;mean(C2)+0.005*rand(1,3)];
    C3 = [C3;mean(C3)+0.005*rand(1,3)];
    CKai = [CKai;mean(CKai)+0.005*rand(1,3)];
    CMir = [CMir;mean(CMir)+0.005*rand(1,3)];
    
    % Get the minimal bounding sphere radii of each algorithm:
    [BST(i,1),~] = ApproxMinBoundSphereND(C1);
    [BST(i,2),~] = ApproxMinBoundSphereND(C2);
    [BST(i,3),~] = ApproxMinBoundSphereND(C3);
    [BST(i,4),~] = ApproxMinBoundSphereND(CKai);
    [BST(i,5),~] = ApproxMinBoundSphereND(CMir);
    
    clearvars C1 C2 C3 CKai CMir
    
    %% PAtella
    
    C1 =  [Results(i, 1).Pat.Center'; Results(i, 2).Pat.Center'; Results(i, 3).Pat.Center'];
    C2 =  [Results(i, 1).Pat.Center3; Results(i, 2).Pat.Center3; Results(i, 3).Pat.Center3];
    C3 =  [Results(i, 1).Pat.Center4; Results(i, 2).Pat.Center4; Results(i, 3).Pat.Center4];
    
    % Because the ApproxMinBoundSphereND function need 4 non-coplanar points
    % a 4th point is added, it is the barycenter of the 3 origin points
    % randomly perturbated in space :
    C1 = [C1;mean(C1)+0.005*rand(1,3)];
    C2 = [C2;mean(C2)+0.005*rand(1,3)];
    C3 = [C3;mean(C3)+0.005*rand(1,3)];
    C0 = C1; % The Center is identically defined in VR and Rainbow 2013
    
    % Get the minimal bounding sphere radii of each algorithm:
    [BSP(i,1),~] = ApproxMinBoundSphereND(C1);
    [BSP(i,2),~] = ApproxMinBoundSphereND(C2);
    [BSP(i,3),~] = ApproxMinBoundSphereND(C3);
    [BSP(i,4),~] = ApproxMinBoundSphereND(C0);
    
    clearvars C1 C2 C3 C0
    
    
end

%% Write results

%% Femur

Methods = [repmat({'A_PCC'}, [N 1]); repmat({'B_PCS'}, [N 1]); repmat({'C_CE'}, [N 1]); repmat({'D_Kai'}, [N 1]); repmat({'E_Miranda'}, [N 1])];
Rb = BSF(:) ;

mBSF = table(Methods,Rb);

writetable(mBSF);

%% Tibia

Methods = [repmat({'C_CASB'}, [N 1]); repmat({'B_ECASE'}, [N 1]);...
    repmat({'A_PIAASL'}, [N 1]) ; repmat({'D_Kai'}, [N 1]); repmat({'E_Miranda'}, [N 1])];

Rb = BST(:) ;

mBST = table(Methods,Rb);
writetable(mBST);

%% Patella

Methods = [repmat({'C_VR'}, [N 1]); repmat({'B_RL'}, [N 1]); repmat({'A_PIAAS'}, [N 1]); repmat({'D_Rainbow'}, [N 1])];

Rb = BSP(:) ;

mBSP = table(Methods,Rb);
writetable(mBSP);

