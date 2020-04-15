function [T] = meanTemp()
%Definimos T como un vector de 3 entradas ya que sabemos que la ecuación
%tendrá 3 soluciones (2 inestables y 1 estable)

%% Define parameters:
eps=0.62;
Q=342;
theta=(5.67e-08);

%% Plot:
%Sabemos que las soluciones están alrededor de:
%T1~230
%T2~265
%T3~285
%Esto nos sirve porque al ser ecuaciones no polinomicas, hay que
%especificarle a Matlab alrededor de donde estan nuestras soluciones; de lo
%contrario, solo nos regresa la primera que encuentre.
syms t
eqnLeft = (1-((.5-(0.2)*(tanh((t-265)/10)))))*Q;
eqnRight = (eps)*(theta)*(t^4);
fplot([eqnLeft eqnRight])
xlim([200 310])
ylim([50 350])
title([texlabel(eqnLeft) ' = ' texlabel(eqnRight)])
xlabel('Temperatura (K)')
ylabel('Energía (Wm^-2)')

%% Solve:

t1 = vpasolve(eqnLeft == eqnRight, t, 230);
t2=vpasolve(eqnLeft == eqnRight, t, 265);
t3=vpasolve(eqnLeft == eqnRight, t, 285);
T=[t1,t2,t3];
%Así tenemos T nuestro vector con tres soluciones
end