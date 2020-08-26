%Rosenbrock
function [r] = frosenbrock(x)
c = 100; 
r = 0;
for i = 1:length(x)/2
    r = r + c * (x(2 * i) - x(2 * i - 1) ^ 2) ^ 2 + (1 - x(2 * i - 1)) ^ 2; 
end