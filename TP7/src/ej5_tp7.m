function ej5_tp7;clc;close all;
    % parametros del sistema con 2 grados de libertad con amortiguamiento
    l1 = 0.9; % longitud [m]
    l2 = 1.4; % longitud [m]

    m=4000; % masa [kg]
    r = 0.64; % radio [m]
    I = m*r^2; % momento de inercia [kg.m^2]

    m1 = m; % masa [kg]
    m2 = I; % masa [kg]

    k = 20000; % rigidez [N/m]
    k1 = k; % rigidez [N/m]
    k2 = k; % rigidez [N/m]

    M = diag([m1,m2]); % matriz de masa
    K = [
          k1+k2  ,  k2*l2-k1*l1  ;
          k2*l2-k1*l1  , k2*l2^2+k1*l1^2
        ]; % matriz de rigidez
    
    c=2000; % amortiguamiento [N.s/m]
    c1 = c; % amortiguamiento [N.s/m]
    c2 = c; % amortiguamiento [N.s/m]

    C = [
          c1+c2  ,  c2*l2-c1*l1  ;
          c2*l2-c1*l1  , c2*l2^2+c1*l1^2
        ]; % matriz de amortiguamiento

    %Condiciones iniciales
    x0 = [0.05 0]'; % desplazamiento inicial [m]
    dx0 = [0 0]'; % velocidad inicial [m/s]

    % Tiempo de analisis
    tf = 10; % tiempo final [s]
    dt = 1e-6; % paso de tiempo [s]
    t = 0:dt:tf; % vector de tiempo [s]

    % autovalores y autovectores
    [V,lambda] = eig(K,M);
    wn =diag(sqrt(lambda)) % frecuencias naturales [rad/s]

    x = Descomposicion_Modal_Amortiguado( M,C,V,wn, x0, dx0, t);

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo con amortiguamiento');
    hold on;
    plot(t,x(1,:),"LineWidth",1.5,"Color","b");
    plot(t,x(2,:),"LineWidth",1.5,"Color","r");
    hold off;
    ylabel("Poscion");
    xlabel("Tiempo [s]");
    legend("Desplazamiento","Rotacion");
    grid on;


end

function x = Descomposicion_Modal_Amortiguado( M,C,V,wn, x0, dx0, t)

    % Paso de coordenadas del espacio fisico al espacio modal
    y0 = V'*M*x0; % desplazamiento inicial en el espacio modal
    dy0 = V'*M*dx0; % velocidad inicial en el espacio modal
    
    C_modal = V'*C*V; % matriz de amortiguamiento modal [kips.sec/in]

    zitta = diag(C_modal)./(2*wn); % relacion de amortiguamiento

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