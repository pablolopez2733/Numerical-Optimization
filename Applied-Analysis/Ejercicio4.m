clear;   close all;   clc;

f = @(X,Y)(1 + ((X + Y + 1).^2) * (19 - (14 * X) + (3 * (X .^2)) - 14*Y + (6 .* X.*Y) + (3 * (Y.^2)))) .* ...
        (30 + ((2 * X - 3 * Y).^2) .* (18 - 32 * X + 12 * (X .^2) + 48 * Y - (36 .* X.*Y) + (27 * (Y.^2))) );
stepsize =  0.1;  % mas chico hace los conjuntos de nivel mas detallados
x0=[1;1];
itmax=500;
delta=1;

[X,Y] = meshgrid(-2:stepsize:2);
z = f(X,Y);
niveles = [0.1, -2:2];
contour(X,Y,z, niveles)
f=@(x)((1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*(30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2)));
axis equal
%Hacemos un arreglo para guardar las iteraciones:
x=zeros(10,1);
y=zeros(10,1);
%iniciamos la primera iteracion en el punto x0
x(1)=x0(1);
y(1)=x0(2);
j=2;
for i=1:6
    [x1, msg] = mRC2(f, x0, itmax);
    x(j)=x1(1);
    y(j)=x1(2);
    j=j+1;
end
hold on
x = x;
y = y;
plot(x,y,'--d')
hold off
