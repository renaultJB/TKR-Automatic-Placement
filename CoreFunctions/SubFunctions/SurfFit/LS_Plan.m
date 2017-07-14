function [Nrml,d] = LS_Plan( Points )
%LS_Plan associate a Least-Square plane on the point cloud Points
%   this is done by computed the covariance matrix of the point cloud and
% 	the eigen vector correspondign to the first eigen value is the normal
%	of the least square plan

[V,~] = eig(cov(Points));

Nrml = V(:,1);
d = -mean(Points)*Nrml;

% Change orientation such that the normal is positive relative to the world coordinate system
[~,Im]=max(abs(Nrml));
if Nrml(Im)<0
    Nrml=-Nrml;
    d=-d;
end


end

