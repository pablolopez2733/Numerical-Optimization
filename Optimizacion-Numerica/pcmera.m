function [x,lambda] = pcmera(Q,A,c,b)
%Método del rango para el problema
%MIN (1/2)x'Qx+c'x
%SA ax = b
%
%Q es nxn spd en Rn y cuya inversa se calcula fácil
%A matriz mxm de rango m.
%c vector columna de orden n
%b vector columna de orden n
%OUT
%x solución del problema
%lamda vector columna de orden m y es el multiplicador de Lagrange.
%--------------------------------------------------------------------------
S = inv(Q);
B = A*S*A';
ld = -(b + A*S*c);

%Usamos GC para resolver
[lambda] = pcg(B,ld);

x = -S * (c+A'*lambda);
