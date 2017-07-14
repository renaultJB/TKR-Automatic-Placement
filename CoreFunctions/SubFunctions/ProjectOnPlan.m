function [ Pts_Proj ] = ProjectOnPlan( Pts , n , Op )
%ProjectOnPlan : project points Pts on a plan describe by its normal n and
%its origin Op (or d the last parameter of the plan equation)
if length(Op)==1
    Op = [0 0 -Op(1)/n(3)]; %Origin points of plan
end
OpPts = bsxfun(@minus,Pts,Op); %Create vector Pts origin of plan Op
Pts_Proj = Pts - (OpPts*n)*n'; % Substract Pts
end

