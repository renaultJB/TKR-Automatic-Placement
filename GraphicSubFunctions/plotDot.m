function plotDot( center, color, r )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin <2
    color='k';
    r=1.75;
end

if nargin <3
    r=1.75;
end

[x,y,z] = sphere(50);
x0 = center(1); y0 = center(2); z0 = center(3);
x = x*r + x0;
y = y*r + y0;
z = z*r + z0;

hold on
% lightGrey = 0.8*[1 1 1]; % It looks better if the lines are lighter
surface(x,y,z,'FaceColor', color,'EdgeColor','none')
hold on

end

