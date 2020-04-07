clear;   close all;   clc;

%% Define f with two arguments  and  so that
% it can be evaluated for matrices of values X, Y


f = @(X,Y)((1+(X+Y+1).^2*(19-14.*X+3.*X.^2-14.*Y+6.*X.*Y+3.*Y.^2))*(30+(2.*X-3.*Y).^2.*(18-32.*X+12.*X.^2+48.*Y-36.*X.*Y+27.*Y.^2)));

%% define point and trust region radius
x0    = [-0.5;-0.5];
delta = 1.5;
pg=@(x)((1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*(30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2)));
g=apGrad(pg,x0);
B=apHess(pg,x0);

%% plot f in cartesian coordinates arround x0
showPlot = true;
if showPlot
	Delta = 1.1*delta;
	uniGrid = linspace(x0(1)-Delta, x0(1)+Delta, 32);
	[X,Y]   = meshgrid(uniGrid, uniGrid);
	Z  = f(X,Y);
	s2 = surf(X,Y,Z);
end


%% plot quadratic model in trust region with polar coordinates arround x0
hold on
% polar coordinates arround x0
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,delta,16));
X     = R.*cos(T) +x0(1);
Y     = R.*sin(T) +x0(2);
% the quadratic model (simple in this case)
%m  = @(X, Y) (1827636-2186063.99970461.*X+6920496.00541966.*Y+(1/2)*((521387.970703125).*X.^2+(2*(-5880852.31494141)).*X.*Y+21866510.1406250.*Y.^2)) ;
m  = @(X, Y) (pg(x0)+(g(1)).*X+(g(2)).*Y+(1/2)*((B(1,1)).*X.^2+((B(1,2))+(B(2,1))).*X.*Y+(B(2,2)).*Y.^2)) ;
% evaluation and plot
Z  = m(X,Y);
s1 = mesh(X,Y,Z);
view(130, 15)
hold off