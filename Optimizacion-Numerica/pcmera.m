function [x,lambda] = pcmera(Q,A,c,b)
%M�todo del rango para el problema
%MIN (1/2)x'Qx+c'x
%SA ax = b
%
%Q es nxn spd en Rn y cuya inversa se calcula f�cil
%A matriz mxm de rango m.
%c vector columna de orden n
%b vector columna de orden n
%OUT
%x soluci�n del problema
%lamda vector columna de orden m y es el multiplicador de Lagrange.
%--------------------------------------------------------------------------
S = inv(Q);
B = A*S*A';
ld = -(b + A*S*c);

%Usamos GC para resolver
[lambda] = pcg(B,ld);

x = -S * (c+A'*lambda);
