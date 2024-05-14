function ej3_tp7;clc;close all;
    % Parametros del sistema con 3 grados de libertad sin amortiguamiento
    m = 4; % Masa [kg]
    m1 = m; % masa [kg]
    m2 = 2*m; % masa [kg]
    m3 = 2*m; % masa [kg]

    k = 4; % Rigidez [N/m]
    k1 = k; % rigidez [N/m]
    k2 = 2*k; % rigidez [N/m]
    k3 = k; % rigidez [N/m]

    M = diag([m1,m2,m3]); % matriz de masa
    K = [k1+k2,  -k2 ,    0;
    -k2   ,  k2+k3 ,  -k3;
    0      , -k3 ,    k3];

    zitta = [0.1 0.1 0.1]; % relacion de amortiguamiento

    % autovalores y autovectores
    [V,lambda] = eig(K,M);

    wn =diag(sqrt(lambda)) % frecuencias naturales [rad/s]
    V_N = V./V(1,:) % autovectores normalizados

    % Condiciones iniciales
    x0 = [1 0 0]'; % desplazamiento inicial [m]
    dx0 = [0 0 0]'; % velocidad inicial [m/s]

    % Tiempo de analisis
    tf = 2; % tiempo final [s]
    dt = 1e-6; % paso de tiempo [s]
    t = 0:dt:tf; % vector de tiempo [s]

    
    x1 = Descomposicion_Modal_Amortiguado(zitta, M,V,wn, x0, dx0, t); % Calculo de la respuesta en el tiempo

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo sin amortiguamiento');
    hold on;
    plot(t,x1(1,:),"LineWidth",1.5,"Color","b");
    plot(t,x1(2,:),"LineWidth",1.5,"Color","r");
    plot(t,x1(3,:),"LineWidth",1.5,"Color","g");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)");
    grid on;


end

function x = Descomposicion_Modal_Amortiguado(zitta, M,V,wn, x0, dx0, t)

    % Paso de coordenadas del espacio fisico al espacio modal
    y0 = V'*M*x0; % desplazamiento inicial en el espacio modal
    dy0 = V'*M*dx0; % velocidad inicial en el espacio modal
    
    % C_modal = 2*zitta.*wn.*eye(size(M)); % matriz de amortiguamiento modal [kips.sec/in]

    wd = wn*sqrt(1-zitta.^2); % Frecuencia de amortiguamiento [rad/s]
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