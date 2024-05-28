function ej3_p;clc;close all
    %Datos
    m=30;
    l=1;
    k=4*10^3;

    %Elementos masa
    m1=m*l^2/9;
    m2=m;
    m3=2*m;

    M=diag([m1,m2,m3]); %Matriz masa

    %Matriz rigidez
    K=[10*k*l^2/9   k*l/3    -2*k*l/3;
    k*l/3        k           0;
    -2*k*l/3     0           k];

    %Tiempo dato
    tf=6;
    dt=0.02;
    t=0:dt:tf;

    %%Condiciones Iniciales

    angulo0_deg=30; %dato
    angulo0_rad = angulo0_deg*pi/180; %angulo en radianes

    xo=[angulo0_rad;0;0]; %dato
    dxo=[0;0;0]; %Velocidades inciales

    %Relacion de amortiguamiento
    zitta=[0 0 0];

    %Autovalores y autovectores
    [V,lambda]=eig(K,M);
    wn=diag(sqrt(lambda));
    V_N=V./V(1,:);

    X1=Descomposicion_Modal_Amortiguado(zitta,M,V,wn,xo,dxo,t);

    %Mostramos la respuesta en el tiempo
    figure('Name','Respuesta en el tiempo sin amortiguamiento');
    hold on;
    plot(t,X1(1,:),'LineWidth',1.5,'Color','b');
    plot(t,X1(2,:),'LineWidth',1.5,'Color','r');
    plot(t,X1(3,:),'LineWidth',1.5,'Color','g');
    hold off;
    ylabel('Desplazamiento [rad]');
    xlabel('Tiempo [s]');
    legend('x1(t)','x2(t)','x3(t)');
    grid on;

    %Calculo de amplitud máxima de theta entre los 0.5 y 1.5 segundos 

indice1=find(t==0.5);
indice2=find(t==1.5);

t_acotado= t(indice1:indice2);

X1_acotado=X1(:,indice1: indice2);

disp("amplitud positiva (máx. valor) de movimiento del grado de libertad, θ(t), dentro del intervalo de tiempo [0.5, 1.5]s. ")
X1_ang_deg = X1_acotado(1,:) .* 180/pi; %a los ángulos los convierto de rad a degree
[theta_max,indice_max] = max(X1_ang_deg)
t_max=t_acotado(indice_max)

%Entre los 2.0 y 3.0 segundos
indice1=find(t==2);
indice2=find(t==3);

disp("amplitud positiva (máx. valor) de movimiento del grado de libertad, x1(t), dentro del intervalo de tiempo [2.0, 3.0]s. ")
t_acotado= t(indice1:indice2);
X1_acotado=X1(:,indice1: indice2);
[x1_max,indice_max] = max(X1_acotado(2,:))
t_max=t_acotado(indice_max)

disp("amplitud positiva (máx. valor) de movimiento del grado de libertad, x2(t), dentro del intervalo de tiempo [2.0, 3.0]s. ")
t_acotado= t(indice1:indice2);
X1_acotado=X1(:,indice1: indice2);
[x2_max,indice_max] = max(X1_acotado(3,:))
t_max=t_acotado(indice_max)


end

function x = Descomposicion_Modal_Amortiguado(zitta, M,V,wn, x0, dx0, t)

    % Paso de coordenadas del espacio fisico al espacio modal
    y0 = V'*M*x0; % desplazamiento inicial en el espacio modal
    dy0 = V'*M*dx0; % velocidad inicial en el espacio modal
    
    % C_modal = 2*zitta.*wn.*eye(size(M)); % matriz de amortiguamiento modal [kips.sec/in]

    wd = wn.*sqrt(1-zitta.^2); % Frecuencia de amortiguamiento [rad/s]
    e = exp(1); 

    Y = zeros(length(wn),length(t)); % desplazamiento en el espacio modal

    for i = 1:length(wn)

        % Primer termino de la solucion homogenea
        Y(i,:) = e.^(-zitta(i)*wn(i)*t).*(cos(wd(i)*t)+zitta(i)/sqrt(1-zitta(i)^2)*sin(wd(i)*t)).*y0(i);
        % Segundo termino de la solucion homogenea
        Y(i,:) = Y(i,:) + e.^(-zitta(i)*wn(i)*t).*(sin(wd(i)*t)).*(dy0(i)/(wd(i)));
    end

    % Paso de coordenadas del espacio modal al espacio fisico
    x = V*Y; % desplazamiento en el espacio fisico

end