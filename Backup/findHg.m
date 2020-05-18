function [p] = findHg(S,Gamma,gnew)
p = gnew;
rhoinv = dot(S,Gamma,1);
alpha = zeros(length(S(1,:)),1);
for i = 1:length(S(1,:))
    alpha(i) = (dot(S(:,i),p)/rhoinv(i));
    p = p - Gamma(:,i)*alpha(i);
end
D = rhoinv(1)/dot(Gamma(:,i),Gamma(:,i));
p = D * p;

for i=1:length(S(1,:))
    B = dot(Gamma(:,i),p)/rhoinv(i);
    p = p + (alpha(i) - B)*S(:,i);
end
end

