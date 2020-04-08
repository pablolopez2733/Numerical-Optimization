%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ejercicio 2
%Dibujar direcciones en 2d. Sea f : R2 a R y el punto inicial x0 cerca 
%de un minimo local x0. La matriz Hessiana en x0 sea simetrica y 
%(semi-)definida positiva y se usa para definir el modelo cuadratico en x0. 
%Escoge una region de confianza con delta > 0. 
%Luego haga un plot (en dos dimensiones) que contiene:
%la frontera de la region de confianza, 
%algunos conjuntos de nivel en R2 del modelo cuadratico en la region de la confianza. 
%los tres direcciones Newton, Cauchy, dogleg. Para obtenerlas use sus 
%funciones pDogLeg, pCauchy.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Limpiamos 
close all
clc
format long
%% 1. Definimos nuestras variables a utilizar
delta = 1;
x0 = [-0.3;-0.7];
f=@(x) (1 + ((x(1) + x(2) + 1).^2) * (19 - (14 * x(1)) + (3 * (x(1) .^2)) - 14*x(2) + (6 .* x(1).*x(2)) + (3 * (x(2).^2)))) .* ...
        (30 + ((2 * x(1) - 3 * x(2)).^2) .* (18 - 32 * x(1) + 12 * (x(1) .^2) + 48 * x(2) - (36 .* x(1).*x(2)) + (27 * (x(2).^2))));
gk = apGrad(f,x0);
Bk = apHess(f,x0);
%% 2. Definimos puntos notables
pN = (-1)*inv(Bk)*gk;
pC = pCauchy(Bk,gk,delta);
pDL = pDogLeg(Bk,gk,delta);
%% 3.Utilizamos el codigo scrLevel set en el modelo cuadrático
mk = @(x,y) f(x0)+(gk(1)).*x+(gk(2)).*y+(1/2)*((Bk(1,1)).*(x.^2)+((Bk(1,2))+(Bk(2,1))).*x.*y+(Bk(2,2)).*(y.^2));
stepsize =  0.01;  
[X,Y] = meshgrid(-2:stepsize:2);
Z = mk(X,Y);
niveles = [-5,0.01:5];
contour(X,Y,Z,niveles);

%% 4. Necesitamos definir una circunferencia
T = linspace(0, 2*pi, 128);
xpolar = delta*cos(T) + x0(1);
ypolar = delta*sin(T) + x0(2);
axis equal
%% 5. Graficamos
hold on

dirN = quiver(x0(1), x0(2),pN(1), pN(2),'Color', 'red',   'LineWidth', 2.5);

dirC = quiver(x0(1), x0(2), pC(1), pC(2),'Color', 'blue',   'LineWidth', 2.5);

dirDog = quiver(x0(1), x0(2), pDL(1), pDL(2),'Color', 'green', 'LineWidth', 2.5);

regionConfianza = plot(xpolar, ypolar,'r--', 'LineWidth',1);

legend([dirN, dirC, dirDog, regionConfianza], {'Direccion Newton', ... 
    'Direccion Cauchy', 'Direccion Dogleg', 'Region de Confianza'});

hold off
grid on
view(2);







