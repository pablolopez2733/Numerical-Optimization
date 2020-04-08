clear;   close all;   clc;

%% Define f with two arguments  and  so that
% it can be evaluated for matrices of values X, Y


f = @(X,Y)(1 + ((X + Y + 1).^2) * (19 - (14 * X) + (3 * (X .^2)) - 14*Y + (6 .* X.*Y) + (3 * (Y.^2)))) .* ...
        (30 + ((2 * X - 3 * Y).^2) .* (18 - 32 * X + 12 * (X .^2) + 48 * Y - (36 .* X.*Y) + (27 * (Y.^2))) );

%% define point and trust region radius
x0    = [0.5;-0.5];
delta = 3.5;
%Now we define our Price-Goldstein function:
f1=@(x) (1 + ((x(1) + x(2) + 1).^2) * (19 - (14 * x(1)) + (3 * (x(1) .^2)) - 14*x(2) + (6 .* x(1).*x(2)) + (3 * (x(2).^2)))) .* ...
        (30 + ((2 * x(1) - 3 * x(2)).^2) .* (18 - 32 * x(1) + 12 * (x(1) .^2) + 48 * x(2) - (36 .* x(1).*x(2)) + (27 * (x(2).^2))));
%We approxiate the gradient and Hessian matrix at x0
%g=apGrad(f1,x0);
%B=apHess(f1,x0);

% plot f in cartesian coordinates arround x0
showPlot = true;
hold on
if showPlot
	Delta = 1.1*delta;
	uniGrid = linspace(x0(1)-1.2*Delta, x0(1)+1.2*(Delta), 32);
	[X,Y]   = meshgrid(uniGrid, uniGrid);
	Z  = f(X,Y);
	s2 = surf(X,Y,Z);
    %zlim([-inf 1000000]);
end


%% plot quadratic model in trust region with polar coordinates arround x0
hold on
% polar coordinates arround x0
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,delta,16));
X = R.*cos(T)+x0(1);
Y = R.*sin(T)+x0(2);
% the quadratic model :
%m  = @(X,Y) f1(x0) +15.*X -22.5.*Y + 0.5.*(B(1,1).*X.^2)+0.5.*((B(1,2)+B(2,1)).*X.*Y)+0.5.*(B(2,2).*Y.^2) ;
m  = @(X,Y) (1.937500000000000e+02) +(-0667.499999692550).*X +(1.582499999511301e+03).*Y + 0.5.*((-5.962497473822700e+02).*X.^2)+0.5.*(((-2.096250659094917e+03)+(-2.096250659094917e+03)).*X.*Y)+0.5.*((6.903749645657010e+03).*Y.^2) ;
% evaluation and plot
Z  = m(X,Y);
s1 = mesh(X,Y,Z);
%zlim([-inf 10000]);
view(130, 15);
hold off
