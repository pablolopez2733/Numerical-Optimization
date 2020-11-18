function[x] = metespnulo(Q,A,c,b)
%------------------------------------------
%Método del espacio nulo
%MIN (1/2)x'Qx+c'x
%SA Ax = b
%
%Q es spd en Rn
%A matriz mxm de rango m.
%c vector columna de orden n
%b vector columna de orden ,
%OUT
%x aprox al mínimo del problema
%-------------------------------------------
Z = null(A);
xp = A\b;
%cambio de variable x = Zy + xp
B = Z'*Q*Z;
d = Z'*(Q*xp+c);
%sistema lineal = By = d
y = -B\d;
x = Z*y + xp;

end