disp(eye(4))

x=[0:0.1:10];

x2=linspace(0,10,5);

result = myFun(x);





% Definir la función de la ecuación diferencial
dydt = @(t,y) y - t^2 + 1;

% Definir el rango de tiempo para la solución
tspan = [0 2];

% Definir la condición inicial
y0 = 0.5;

% Llamar a ode45 para resolver la ecuación diferencial
[t,y] = ode45(dydt, tspan, y0);

% Graficar la solución
plot(t,y)



function y = myFun(x)

    y = x.^2;

end