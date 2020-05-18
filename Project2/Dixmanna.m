function [f] = Dixmanna(x0)
%% Parameters (L)
alpha=1;
beta=.26;
sigma=.26;
gamma=.26;
k1=2;
k2=0;    
k3=0;
k4=2;

n=length(x0);
m=floor(n/3);
term1=0;
term2=0;
term3=0;
term4=0;
x=x0;

%% Dixmanna function:
%termino1
for i=1:n
    term1= term1 + (alpha)*(x(i)^2)*(i/n)^k1;
end
%termino2
for i=1:(n-1)
    term2= term2 + (beta)*(x(i)^2)*(x(i+1)+x(i+1)^2)^2*(i/n)^k2;
end
%termino3
for i=1:2*m
    term3= term3 + (gamma)*(x(i)^2)*(x(i+m)^4)*(i/n)^k3;
end
%termino4
for i=1:m
    term4= term4 + (sigma)*(x(i))*(x(i+2*m))*(i/n)^k4;
end

f=1+term1+term2+term3+term4;

end