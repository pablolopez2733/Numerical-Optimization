clear;
close all;
clc;

%Graficamos las curvas de nivel
f = @(x,y)  ((1)+(((x+y+1).^2).*(19-14.*x+3.*x.*x-14.*y+6.*x.*y+3.*y.*y))).*((30)+(((2.*x-3.*y).^2).*(18-32.*x+12.*x.*x+48.*y-36.*x.*y+27.*y.*y))) ;
stepsize=0.01;
[X,Y]=meshgrid(-2.7:stepsize:2.7);
z=f(X,Y);

contour(X,Y,z,40)

axis equal

%Graficamos el recorrido de las iteraciones
hold on
x0=[0;-2.6];
XY=zeros(25,2);
XY(1,:)=x0;
j=2;
for i=1:24
    [x,msg,it]= mRC1(@(x) ((1)+(((x(1,1)+x(2,1)+1)^2)*((19)+(-14*x(1,1))+(3*x(1,1)*x(1,1))+(-14*x(2,1))+(6*x(1,1)*x(2,1))+(3*x(2,1)*x(2,1)))))*((30)+((((2*x(1,1))+(-3*x(2,1)))^2)*((18)+(-32*x(1,1))+(12*x(1,1)*x(1,1))+(48*x(2,1))+(-36*x(1,1)*x(2,1))+(27*x(2,1)*x(2,1))))),x0,i);
    XY(j,:)=x;
    j=j+1;
end

plot(XY(:,1),XY(:,2),'--d')
hold off
