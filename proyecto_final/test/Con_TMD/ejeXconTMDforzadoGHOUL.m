function ejeXsinTMDforzadoGHOUL

clc; clear all;

%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
m1=62000; %% Masa 1
m2=93000; %% Masa 2
m3=0.03*(m1+m2); %% Masa TMD

z=[0.015;0.015;0.01];

k1=2e8; %% Rigidez 1
k2=9e6; %% Rigidez 2
k3=(9.4)^2*m3 %% Rigidez 3

M= [m1 0 0;
    0 m2 0;
    0 0 m3]; %% Matriz de masa

K=[k1+k2 -k2 0;
    -k2   k2 -k3;
       0 -k3 k3]; %% Matriz de Rigidez

dt = 0.005;
t = 0:dt:50;
length(t)

xINI = [0; 0; 0];
dxINI = [0; 0; 0];

%% ENCUENTRO AUTOVECTORES Y AUTOVALORES
[X, lambda] = eig(K, M); % el X que me larga ya esta normalizado de modo q X'*m*X=I

wn = diag(sqrt(lambda));
wd = wn .* sqrt(1 - z.^2);

%% MÉTODO DE DESCOMPOSICIÓN MODAL

Mmodal = round(X' * M * X); % M modal
Kmodal = round(X' * K * X); % K modal
Cmodal = 2 * z .* wn .* eye(size(Mmodal)); % C modal

Y0 = X' * M * xINI;
dY0 = X' * M * dxINI; % paso las condiciones inciales a coord modales

% Ahora el sistema será Mmodal*y''+Cmodal*y'+Kmodal*y=0 es un sistema de EDOs
% desacoplado y habra wd
wd = wn .* sqrt(1 - z.^2);

% Transitoria:

Y = zeros(length(wn), length(t));

for i = 1:length(wn)

    % 1er término
    Y(i, :) = exp(-z(i) * wn(i) * t) .* (cos(wd(i) * t) + z(i) / (sqrt(1 - z(i)^2)) .* sin(wd(i) * t)) .* Y0(i);

    % 2do término
    Y(i, :) = Y(i, :) + (1 / wd(i) .* exp(-z(i) * wn(i) * t) .* sin(wd(i) * t)) .* dY0(i);

end

%% PASO DE COORD MODALES A GEOMETRICAS

xt = X * Y;

%% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA

Twa = 5; %% Periodo del oleaje
VAmax = 2; %% Velocidad máxima del agua
Cd = 0.7; %% Coeficiente de arrastre
RoW = 1000; %% Densidad del agua
Ainmersa = 14 * 4.176; %% Área a considerar, h sobre agua 14m y diametro 4.176m

VAt = zeros(1, length(t)); %% Vector agua
for i = 1:length(t)

    tmod = mod(t(i), Twa); % llevo el valor sub i de t a un valor dentro del periodo (entre 0 y Tp)
    VAt(i) = (tmod >= 0) .* (tmod <= Twa * 4 / 5) .* (((VAmax / 2) / (Twa * 4 / 5)) * tmod + VAmax / 2) + (tmod > Twa * 4 / 5) .* (tmod <= Twa) .* (-VAmax / 2 * tmod + 6);
end

subplot(3,1,1)
plot(t, VAt, 'LineWidth', 2)
title('Velocidad del agua [m/s]')

FW = (0.5 * RoW * Cd * Ainmersa .* VAt.^2 + 91120);

subplot(3,1,2)
plot(t, FW, 'LineWidth', 2)
title('Fuerza del agua [N]')

%% DEFINO LA VELOCIDAD DEL VIENTO Y PARAMETROS DEL VIENTO


Aviento=2000;
frecViento=0.5;

%Fviento = 100000 * ones(1, length(t));

%Agrego variacion al viento
for i=1:2500
    Fviento(i)=100000+Aviento*sin(2*pi*frecViento*t(i))+Aviento/3*sin(3*2*pi*frecViento*t(i))+Aviento/4*cos(4.1*2*pi*frecViento*t(i))+Aviento/8.2*sin(8.34*2*pi*frecViento*t(i));
end
for i=2501:5000
    Fviento(i)=100000+Aviento/1.2*sin(1.1*2*pi*frecViento*t(i))+Aviento/3.23*sin(4.21*2*pi*frecViento*t(i))+Aviento/5.23*cos(6.23*2*pi*frecViento*t(i))+Aviento/8.26*sin(5.32*2*pi*frecViento*t(i));
end
for i=5001:7500
    Fviento(i)=100000+Aviento/1.77*cos(1.1*2*pi*frecViento*t(i))+Aviento/3.91*sin(2.21*2*pi*frecViento*t(i))+Aviento/3.678*cos(7.1*2*pi*frecViento*t(i))+Aviento/9.4*sin(4.375*2*pi*frecViento*t(i));
end
for i=7501:10001
    Fviento(i)=100000+Aviento/0.992*cos(1.632*2*pi*frecViento*t(i))+Aviento/5.13*sin(4.581*2*pi*frecViento*t(i))+Aviento/6.8*sin(7.1*2*pi*frecViento*t(i))+Aviento/9.4*sin(4.65*2*pi*frecViento*t(i));
end

subplot(3,1,3)
plot(t, Fviento, 'LineWidth', 2)
title('Fuerza del viento [N]')

%% RESPUESTA PERMANENTE

F = zeros(3, length(t)); % Vector de fuerzas
F(1, :) = FW;            % Fuerza periódica en el grado de libertad 1
F(2, :) = Fviento;       % Fuerza constante en el grado de libertad 2

% Transformar fuerzas a coordenadas modales
Fmodal = X' * F;

% Frecuencias de Fourier
omega = 2 * pi * (0:(length(t)-1)) / (length(t) * dt);

% Función de transferencia H(iw)
H = zeros(size(Mmodal, 1), length(omega));
for i = 1:length(omega)
    H(:, i) = diag(((-omega(i)^2 * Mmodal + 1i * omega(i) * Cmodal + Kmodal) \ eye(size(Mmodal))));
end

% Respuesta en frecuencia
Y_freq = H .* fft(Fmodal, [], 2);

% Transformar de vuelta al dominio del tiempo
Y_perm = ifft(Y_freq, [], 2, 'symmetric');

% Transformar respuesta permanente a coordenadas físicas
x_perm = X * Y_perm;

xGeneral=xt+x_perm;

%% Ploteado
figure(2)
plot(t, xGeneral(1, :), 'b', 'LineWidth', 2)
hold on 
plot(t, xGeneral(2, :), 'r', 'LineWidth', 2)
plot(t, xGeneral(3, :), 'k', 'LineWidth', 2)
hold off
grid on
title('Respuesta Permanente')
legend('Permanente GDL 1', 'Permanente GDL 2')

