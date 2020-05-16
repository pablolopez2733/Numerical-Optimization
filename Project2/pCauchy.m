function [pC] = pCauchy( B, g, delta)


pk=(delta/norm(g))* g;



%Check if g*B*g is <=0

tao=0;

if(g.'*B*g <=0)

    tao=1;

else

    tao=min((norm(g)^3)/(delta*(g.'*B*g)),1);

    

end

pC=(-1)*(tao)*pk;





end