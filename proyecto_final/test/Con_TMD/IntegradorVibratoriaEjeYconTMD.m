function IntegradorVibratoriaEjeYconTMD

clc; clear all; close all;

%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
m1=62000; %% Masa 1
m2=93000; %% Masa 2
m3=0.05*(m1+m2); %% Masa TMD

z=[0.015;0.015;0.1];

k1=2e8; %% Rigidez 1
k2=9e6; %% Rigidez 2
k3=(9)^2*m3 %% Rigidez 3

M= [m1 0 0;
    0 m2 0;
    0 0 m3]; %% Matriz de masa

K=[k1+k2 -k2 0;
    -k2   k2 -k3;
       0 -k3 k3]; %% Matriz de Rigidez

dt=0.005;
t=0:dt:50;

xINI=[0.4;0;0];
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


figure (1)
plot(t,xt(1,:),'b');
hold on
plot(t,xt(2,:),'r');
hold off

end




