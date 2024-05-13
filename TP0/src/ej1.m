
% Definimos el rango de x
x=[-30:0.5:30];

% Inicializamos y como un vector de ceros del tamaño de x
y = zeros(size(x));

% Definimos las condiciones para y
y(x < -2) = -2;
y(-2 <= x & x < 3) = x(x >= -2 & x < 3).^2;
y(3 <= x & x <= 10) = 1./x(x >= 3 & x <= 10);
y(15 < x & x < 20) = x(x > 15 & x < 20);
y(x == 22) = 3 - x(x == 22);

% Graficamos la función
plot(x,y)


