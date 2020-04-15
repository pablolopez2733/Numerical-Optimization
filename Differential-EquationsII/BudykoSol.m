function [T,E] = BudykoSol()
%Esta función nos regresa los valores de T y E que son soluciones de la
%ecuación para el modelo Lineal Ajustado de Budyko.

%Definimos T como un vector de 3 entradas ya que sabemos que la ecuación
%tendrá 3 soluciones (2 inestables y 1 estable).



%% Define parameters:
% La bibliografía nos dice que:  best fit with observational data
%for the current climate on the Northern Hemisphere is obtained with A= 203.3 Wm-2
% and B = 2.09 Wm-2deg-1 
A= 203.3;
B = 2.09;
Q=342;
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
eqnRight = A+B*t;
fplot([eqnLeft eqnRight])
xlim([-100 310])
ylim([-100 350])
title([texlabel(eqnLeft) ' = ' texlabel(eqnRight)])
xlabel('Temperatura (°C)')
ylabel('Energía (Wm^-2)')

%% Solve:

t1 = vpasolve(eqnLeft == eqnRight, t, 230);
t2=vpasolve(eqnLeft == eqnRight, t, 265);
t3=vpasolve(eqnLeft == eqnRight, t, 285);
T=[t1,t2,t3];
%Así tenemos T nuestro vector con tres soluciones

%% Get Energy values:
E=zeros(3,1);
for i=1:3
    E(i)=A+B*(T(i));
end

end