function [x, lambda] = pc (Q,A,c,b)
%Metodo para resolver:
%MIN (1/2)x'Qx+c'x
%SA ax = b
%-----------------------------------------------------------------------
m = length(b);
n = length(c);
K = [Q A'; A zeros(m)]; 
ld = [-c;b];

%Resolver el sistema
% w = [x;lamda];

w = K\ld;
x = w(1:n);
lambda =  w(n+1:n+m);
end