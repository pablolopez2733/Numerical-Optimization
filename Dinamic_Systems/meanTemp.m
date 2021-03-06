function [T,E] = meanTemp()
%Esta funci�n nos regresa los valores de T y E que son soluciones de la
%ecuaci�n.

%Definimos T como un vector de 3 entradas ya que sabemos que la ecuaci�n
%tendr� 3 soluciones (2 inestables y 1 estable)

%% Define parameters:
eps=0.62;
Q=342;
sigma=(5.67e-08);

%% Plot:
%Sabemos que las soluciones est�n alrededor de:
%T1~230
%T2~265
%T3~285
%Esto nos sirve porque al ser ecuaciones no polinomicas, hay que
%especificarle a Matlab alrededor de donde estan nuestras soluciones; de lo
%contrario, solo nos regresa la primera que encuentre.
syms t
eqnLeft = (1-((.5-(0.2)*(tanh((t-265)/10)))))*Q;
eqnRight = (eps)*(sigma)*(t^4);
fplot([eqnLeft eqnRight])
xlim([200 310])
ylim([50 350])
title([texlabel(eqnLeft) ' = ' texlabel(eqnRight)])
xlabel('Temperatura (K)')
ylabel('Energ�a (Wm^-2)')

%% Solve:

t1 = vpasolve(eqnLeft == eqnRight, t, 230);
t2=vpasolve(eqnLeft == eqnRight, t, 265);
t3=vpasolve(eqnLeft == eqnRight, t, 285);
T=[t1,t2,t3];
%As� tenemos T nuestro vector con tres soluciones

%% Get Energy values:
E=zeros(3,1);
for i=1:3
    E(i)=(T(i)^4)*eps*sigma;
end

end