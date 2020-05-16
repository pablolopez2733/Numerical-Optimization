function [ySol] = TempEq()
eps=0.62;
Q=342;
sigma=(5.67e-08);
syms y(t);
ode = diff(y(t),t) ==((1-((.5-(0.2)*(tanh((y(t)-265)/10)))))*Q)-((eps)*(sigma)*(y(t))^4) ;
ySol(t) = dsolve(ode);

%Try and plot:
fimplicit(@(y,t)ySol,[100,350]);
xlim([200 310])
ylim([50 350])
end