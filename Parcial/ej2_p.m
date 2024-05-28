function ej2_p;clc;close all;
    % Parametros del sistema 
    k1 = 1000; % rigidez [N/m]
    k2 = 500; % rigidez [N/m]
    c = 500; % amortiguamiento [N.s/m]

    m = 10; % masa [kg]
    r = 0.05; % radio [m]
    Jo = 1; % momento de inercia [kg.m^2]

    Fo = 50; % fuerza [N]
    wo = 20; % frecuencia [rad/s]

    dt=0.005; % paso de tiempo [s]
    tf=10; % tiempo final [s]
    t=0:dt:tf; % vector de tiempo [s]

    % Matrices de masa, rigidez y amortiguamiento
    M = diag([3.3333,30,60]);

    K =[
        4.4444e3 1.33333e3 -2.66667e3;
        1.33333e3 4e3 0;
        -2.66667e3 0 4e3];

    % autovalores y autovectores
    [V,lambda] = eig(K,M);
    wn =diag(sqrt(lambda)) % frecuencias naturales [rad/s]

    % Condiciones iniciales
    x0 = [0 0 0]'; % desplazamiento inicial [m]
    dx0 = [0 0 0]'; % velocidad inicial [m/s]

    %

    

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