function IntegradorVibratoriaEjeXconTMDforzado

clc; clear all; close all;

%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
m1=62000; %% Masa 1
m2=93000; %% Masa 2
m3=0.05*(m1+m2); %% Masa TMD

z=[0.015;0.015;0.1];

k1=2e8; %% Rigidez 1
k2=9e6; %% Rigidez 2
k3=(0.5*2*pi)^2*m3 %% Rigidez 3

M= [m1 0 0;
    0 m2 0;
    0 0 m3]; %% Matriz de masa

K=[k1+k2 -k2 0;
    -k2   k2 -k3;
       0 -k3 k3]; %% Matriz de Rigidez

dt=0.005;
t=0:dt:100;

xINI=[0;0;0];
dxINI=[0;0;0];

%% ENCUENTRO AUTOVECTORES Y AUTOVALORES
[X,lambda]=eig(K,M); %el X que me larga ya esta normalizado de modo q X'*m*X=I

wn=diag(sqrt(lambda))
wd=wn.*sqrt(1-z.^2);

%% MÉTODO DE DESCOMPOSICIÓN MODAL

Mmodal=round(X'*M*X); %M modal
Kmodal=round(X'*K*X); %K modal

Y0=X'*M*xINI;
dY0=X'*M*dxINI; %paso las condiciones inciales a coord modales

%Ahora el sistema será Mmodal*y''+Cmodal*y'+Kmodal*y=0 es un sistema de EDOs
%desacoplado y habra wd
wd=wn.*sqrt(1-z.^2)

%Transitoria:

Y=zeros(length(wn),length(t));

for i=1:length(wn)

% 1er término
    Y(i,:)=exp(-z(i)*wn(i)*t).*(cos(wd(i)*t)+z(i)/(sqrt(1-z(i)^2)).*sin(wd(i)*t)).*Y0(i);

    % 2do término
    Y(i,:)=Y(i,:)+(1/wd(i).*exp(-z(i)*wn(i)*t).*sin(wd(i)*t)).*dY0(i);

end

%% PASO DE COORD MODALES A GEOMETRICAS

xt=X*Y;

%% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA

Twa=5; %% Periodo del oleaje
VAmax=2; %% Velocidad máxima del agua
Cd=0.7; %% Coeficiente de arrastre
RoW=1000; %% Densidad del agua
Ainmersa=4*4.176; %% Área a considerar, h sobre agua 14m y diametro 4.176m

VAt=zeros(1,length(t)); %% Vector agua
for i=1:length(t)

    tmod=mod(t(i),Twa); %llevo el valor sub i de t a un valor dentro del periodo (entre 0 y Tp)
    VAt(i)= (tmod>=0).*(tmod<=Twa*4/5).*(((VAmax/2)/(Twa*4/5))*tmod+VAmax/2)+(tmod>Twa*4/5).*(tmod<=Twa).*(-VAmax/2*tmod+6);
end


subplot(3,1,1)
plot(t,VAt,'LineWidth',2)
title('Velocidad del agua [m/s]')

FW=(0.5*RoW*Cd*Ainmersa.*VAt.^2+91120);

subplot(3,1,2)
plot(t,FW,'LineWidth',2)
title('Fuerza del agua [N]')

%% DEFINO LA VELOCIDAD DEL VIENTO Y PARAMETROS DEL VIENTO

frecviento=0.5;
Aviento=1000;

for i=1:length(t)
    Fviento(i,1)=(100000+Aviento*sin(2*pi*frecviento*t(i)));
end

subplot(3,1,3)
plot(t,Fviento,'LineWidth',2)
title('Fuerza del viento [N]')

%% RESPUESTA PERMANENTE

F = zeros(3,length(t)); % Vector de fuerzas
F(1,:) = FW;            % Fuerza del agua en el grado de libertad 2
F(2,:) = Fviento;         % Fuerza de desbalance en el grado de libertad 3

% Transformar fuerzas a coordenadas modales
Fmodal = X' * F;

% Respuesta permanente usando Integral de Duhamel
Y_perm = zeros(length(wn), length(t));
for i = 1:length(wn)
    for j = 1:length(t)
        integral_value = 0;
        for k = 1:j
            h = (1/Mmodal(i,i)) * exp(-z(i) * wn(i) * (t(j) - t(k))) * sin(wd(i) * (t(j) - t(k)));
            integral_value = integral_value + h * Fmodal(i,k) * dt;
        end
        Y_perm(i,j) = integral_value;
    end
end

% Transformar respuesta permanente a coordenadas físicas
x_perm = X * Y_perm;

%% TMD



%% Ploteado
figure(2)
plot(t,x_perm(1,:),'b','LineWidth',2)
hold on 
plot(t,x_perm(2,:),'r','LineWidth',2)
plot(t,x_perm(3,:),'LineWidth',2)
hold off
title('Respuesta Transitoria y Permanente')
legend('Permanente GDL 1','Permanente GDL 2','Permanente GDL 3')

%xgen=xt+x_perm;

% % Cálculo de la FFT de la carga FW
% n = length(FW); % Número de puntos de la señal
% Fs = 1/dt; % Frecuencia de muestreo
% f = (0:n-1)*(Fs/n); % Vector de frecuencias
% 
% FW_fft = fft(x_perm(2,:)); % Transformada de Fourier
% 
% % Solo se considera la mitad del espectro (frecuencias positivas)
% P2 = abs(FW_fft/n);
% P1 = P2(1:floor(n/2)+1); % Asegurando que los índices sean enteros
% P1(2:end-1) = 2*P1(2:end-1);
% f = f(1:floor(n/2)+1); % Asegurando que los índices sean enteros
% 
% % Graficar el espectro de frecuencias
% figure (5)
% plot(f, P1,'LineWidth',2)
% title('Espectro de Frecuencias del gl2')
% xlabel('Frecuencia (Hz)')
% ylabel('|P1(f)|')
% grid on

end




