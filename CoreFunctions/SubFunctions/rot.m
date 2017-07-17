function R = rot(k,fi)

fi = double(fi * pi / 180) ;
k = double(k) / norm(double(k)) ;

% This is just to make it easier to read!
x = k(1);
y = k(2);
z = k(3);

% Create a 3x3 zero matrix
R = zeros(3,3);
% We use the formula for rotating a matrix about a unit vector k

R(1,1) = cos(fi)+x^2*(1-cos(fi));
R(1,2) = x*y*(1-cos(fi))-z*sin(fi);
R(1,3) = x*z*(1-cos(fi))+y*sin(fi);

R(2,1) = y*x*(1-cos(fi))+z*sin(fi);
R(2,2) = cos(fi)+y^2*(1-cos(fi));
R(2,3) = y*z*(1-cos(fi))-x*sin(fi);

R(3,1) = z*x*(1-cos(fi))-y*sin(fi);
R(3,2) = z*y*(1-cos(fi))+x*sin(fi);
R(3,3) = cos(fi)+z^2*(1-cos(fi));
end