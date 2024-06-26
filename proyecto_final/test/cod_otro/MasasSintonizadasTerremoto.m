% PROYECTO FINAL VIBRATORIA MASAS SINTONIZADAS
clc; clear; 

load('VectorAceleraciones.mat', 'a');

% Definicion de variables
m = 36000;
n = 3;
E = 23500 * 10^6;
I = 1.56 * 10^(-4);
h = 3;
k = 2 * 12 * E * I / h^3;
% zita = [0.05; 0.05; 0.05; 0.01; 0.01; 0.01];    %den hartog
g = 9.81;
% F = m * a;
% F = repmat(F, 1, 6)';
% for i = 1:3
%     F(i+3 ,:) = 0;
% end
% Valores iniciales
x0 = [0; 0; 0; 0; 0; 0];
v0 = [0; 0; 0; 0; 0; 0];
% Tiempos
dt = 0.02;
tf = 40;%(length(F)-1)*dt;
t = 0:dt:tf;
g = false;
%% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA
FW=agua(t,g);
    
%% DEFINO LA VELOCIDAD DEL VIENTO Y PARAMETROS DEL VIENTO
% vel_viento=    ; %% Velocidad del viento [nudos]
Fv = viento(t,g);

%% DEFINO PARAMETROS Y FUERZA DE DESBALANCE

Fd = desbalance(t,g);

%% DEFINO LA FUERZA EXTERNA
F = fuerza_externa(t,x0,FW,Fv,Fd);



% Definicion de matrices de masas y rigidez
M_masas = diag(repmat(m, 1, n));
K_rigidez = [2 * k, -k, 0;
            -k, 2 * k, -k;
             0,   -k,   k];

% Calculo de autovalores y autovectores normalizados
[V1, lambda1] = eig(K_rigidez, M_masas);
w1 = sqrt(diag(lambda1));

% Calculo de frecuencias y longitudes de péndulo
f = w1 ./(2 * pi);
l = g ./ (w1.^2);
l(1) = g/(3.5*2*pi)^2;

% Calculo de constantes de los péndulos sintonizados
mp = 0.05 * m * 3;

% Definicion de matrices de masa y rigidez para el sistema sintonizado
% Mp_masas = zeros(size(6, 6));
% for i = 1:3
%     Mp_masas(i, i) = m;
% end
% Mp_masas(6, 6) = mp;
% Mp_masas(5, 5) = 0.1*Mp_masas(6, 6);
% Mp_masas(4, 4) = 0.1*Mp_masas(5, 5);

% kp = zeros(size(3, 1));
% for i = 1:3
%     kp(i, 1) = Mp_masas(i+3, i+3)*g/l(4-i);
% end

% Kp_rigidez = [ 2*k + kp(1), -k, 0, -kp(1), 0, 0;
%               -k, 2*k + kp(2), -k, 0, -kp(2), 0;
%                0, -k, k + kp(3), 0, 0, -kp(3);
%               -kp(1), 0, 0, kp(1), 0, 0;
%                0, -kp(2), 0, 0, kp(2), 0;
%                0, 0, -kp(3), 0, 0, kp(3)];

%% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
m1=62000; %% Masa 1
m2=93000; %% Masa 2

zita=[0.015;0.015;0.015;0.015;0.01;0.01];

k1=2.01e8; %% Rigidez 1
k2=9.36e6; %% Rigidez 2

wdamp=9.7998;
mdamp=0.05*(m1+m2);
k3=(wdamp^2)/mdamp;


Mp_masas=transpose([m1,m1,m2,m2,mdamp,mdamp]).*eye(6);

Kp_rigidez=[k1+k2 0 -k2 0 0 0;
    0 k1+k2 0 -k2 0 0;
    -k2 0 k2+k3 0 -k3 0;
     0 -k2  0 k2+k3 0 -k3;
     0 0 -k3 0 k3 0;
     0 0 0 -k3 0 k3];%% Matriz de Rigidez



% Calculo de autovalores y autovectores normalizados para el sistema sintonizado
[V, lambda] = eig(Kp_rigidez, Mp_masas);
w = sqrt(diag(lambda));



% Normalizar los autovectores
% V_norm = V ./ V(1,:);

% Calculo de masas y rigideces modales
M_modal = V' * Mp_masas * V;
K_modal = V' * Kp_rigidez * V;

% Paso a coordenadas modales las condiciones iniciales
y0 = V' * Mp_masas * x0;
vy0 = V' * Mp_masas * v0;

% CALCULO DE EDO EN VIBRACIONES LIBRES CON PENDULOS CON AMORTIGUACION



% Frecuencias naturales amortiguadas
wd = w .* sqrt(1 - zita.^2);

% Calcular desplazamientos modales en el tiempo
F_modal = V' * F;

yc = zeros(length(w), length(t));
ys = zeros(length(w), length(t));
A = zeros(length(w), length(t));
B = zeros(length(w), length(t));
X = zeros(length(w), length(t));
yt = zeros(length(w), length(t));

for i = 1:length(t)
    for n = 1:length(w)
        yc(n, i) = F_modal(n, i)*cos(wd(n)*t(i));
        ys(n, i) = F_modal(n, i)*sin(wd(n)*t(i));
        if i == 1
            A(n, i) = 0;
            B(n, i) = 0;
        else 
            A(n, i) = A(n, i-1)*exp(-zita(n)*w(n)*dt)+dt/(2*M_modal(n, n)*wd(n))*(yc(n, i-1)*exp(-zita(n)*w(n)*dt)+yc(n, i));
            B(n, i) = B(n, i-1)*exp(-zita(n)*w(n)*dt)+dt/(2*M_modal(n, n)*wd(n))*(ys(n, i-1)*exp(-zita(n)*w(n)*dt)+ys(n, i));
        end
        X(n, i) = A(n, i)*sin(wd(n)*t(i)) - B(n, i)*cos(wd(n)*t(i));
        yt(n, i) = exp(-zita(n) * w(n) * t(i)) * ...
                (cos(wd(n)* t(i)) + zita(n)/sqrt(1 - zita(n)^2)*sin(wd(n)*t(i)))*y0(n) + ...
                (1/wd(n)*exp(-zita(n)*w(n)* t(i))*sin(wd(n)*t(i)))*vy0(n) + X(n, i);
    end
end

% Paso de coordenadas modales a geometricas
xt = V * yt;
respuesta_tmd(t,xt,"RESPUESTA DEL SISTEMA");
% Graficar la solución para vibraciones libres con amortiguación
% figure;
% plot(t, xt(1,:), 'r', t, xt(2,:), 'g', t, xt(3,:), 'b');
% xlabel('Tiempo');
% ylabel('Desplazamientos (ROJO 1, VERDE 2, AZUL 3)');
% title('Respuesta del sistema ');
% grid on;

function FW=agua(t,graph) %% Función que devuelve la velocidad y la fuerza del agua en función del tiempo	
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

    VAt=VAt+VAmax/2;
    FW=0.5*RoW*Cd*Ainmersa.*VAt.^2.*(1-exp(-t./5));

    if graph
        figure(1);
        subplot(4,1,1)
        plot(t,VAt,'LineWidth',2)
        title('Velocidad del agua [m/s]')

        subplot(4,1,2)
        plot(t,FW,'LineWidth',2)
        title('Fuerza del agua [N]')
    end
end

function Fviento=viento(t,graph) %% Función que devuelve la fuerza del viento en función del tiempo
    
    % Agregar el calculo de la fuerza del viento en función del viento
    frecviento=1;
    Aviento=2000;
    for i=1:length(t)
       Fviento(i,1)=(100000+Aviento*sin(2*pi*frecviento*t(i))).*(1-exp(-t(i)/5));
    end

    if graph
        figure(1);
        subplot(4,1,3)
        plot(t,Fviento,'LineWidth',2)
        title('Fuerza del viento [N]')
    end
end

function Fdesb=desbalance(t,graph) %% Función que devuelve la fuerza de desbalance en función del tiempo

    F0Desbalance=1028.93*10;
    wRotor=9.7998;

    Fdesb=F0Desbalance.*sin(wRotor*t);

    if graph
        figure(1);
        subplot(4,1,4)
        plot(t,Fdesb,'LineWidth',2)
        title('Fuerza de desbalance [N]')
    end

end

function P=fuerza_externa(t,x0,Fa,Fv,Fd) %% Función que devuelve la fuerza externa en función del tiempo
    fv_ar = 100000*sin(8.3180*t);

    P=zeros(length(x0),length(t)); % fuerza externa en el tiempo
    P(1,:) = 0; 
    P(2,:) = fv_ar;
    P(3,:) = Fd;
    % P(4,:) = Fv;


end

function respuesta_tmd(t,X,titulo)
    % Mostramos la respuesta en el tiempo
    figure('Name', titulo);
    hold on;
    plot(t,X(1,:),"b");
    plot(t,X(2,:),"r");
    plot(t,X(3,:),"g");
    plot(t,X(4,:),"k");
    plot(t,X(5,:),"m--");
    plot(t,X(6,:),"c--");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)","x4(t), tmdx(t), tmdy(t)");
    grid on;
    
end