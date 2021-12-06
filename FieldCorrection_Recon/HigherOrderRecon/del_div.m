function out = del_div(f)
% to calculate del (div(f)/|div(f)|)
% f is a 2D matrix with the size of nx x ny
[Fx, Fy] =  gradient(f);
absF = sqrt((abs(Fx)).^2+(abs(Fy)).^2);
Fx = Fx./absF;
Fy = Fy./absF;
[Fxx, ~] = gradient(Fx);
[~, Fyy] = gradient(Fy);
out = Fxx+Fyy;
