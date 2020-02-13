function plotCurvature2D( X,C,h )
%PLOTCURVATURE2D plot curvature of the X curve with color C on optional
%figure h
%
% An assumption is made that the curve is described by the 1st 2 columns or
% lines of X, so the code does not work if the curve is in oXZ or oYZ
% planes.
%
%   Inputs:
%   - X a matrix of the 2D curves
%   - C color vector same length as X
%--------------------------------------------------------------------------
if nargin <= 1
    error('Not enough input arguments')
end


% Transpose inputs if necessary
if size(X,1)>size(X,2)
    X = X';
end

if size(C,1)==length(C)
    C = C';
end

x = X(1,:);
y = X(2,:);
z = zeros(size(x));
if nargin==2
    figure()
    surface([x;x],[y;y],[z;z],[C;C],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2);
    axis equal
            
elseif nargin ==3
    h.Visible = 'Off';
    h.Visible = 'On';
    hold on
    surface([x;x],[y;y],[z;z],[C;C],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2);
    axis equal
end


end

