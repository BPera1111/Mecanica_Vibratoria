function ej1_tp7;clc;close all;
    % Parametros del sistema con 3 grados de libertad sin amortiguamiento
    m1 = 2; % masa [kips.sec^2/in]
    m2 = 2; % masa [kips.sec^2/in]
    m3 = 2; % masa [kips.sec^2/in]

    k1 = 600; % rigidez [kips/in]
    k2 = 1200; % rigidez [kips/in]
    k3 = 2400; % rigidez [kips/in]

    M = diag([m1,m2,m3]); % matriz de masa
    K = [
          k1  ,  -k1  , 0    ;
         -k1  , k1+k2 , -k2  ;
          0   ,  -k2  , k2+k3;
        ]; % matriz de rigidez


    % Tiempo de analisis
    tf = 2; % tiempo final [s]
    dt = 1e-6; % paso de tiempo [s]
    t = 0:dt:tf; % vector de tiempo [s]

    % Condiciones iniciales
    x0 = [ 0.3; -0.8 ; 0.3 ]; % desplazamiento inicial [in]
    dx0 = [ 0; 0; 0 ]; % velocidad inicial [in/s]

    x1 = Descomposicion_Modal_No_Amortiguado(M, K, x0, dx0, t); % Calculo de la respuesta en el tiempo

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo sin amortiguamiento');
    hold on;
    plot(t,x1(1,:),"LineWidth",1.5,"Color","b");
    plot(t,x1(2,:),"LineWidth",1.5,"Color","r");
    plot(t,x1(3,:),"LineWidth",1.5,"Color","g");
    hold off;
    ylabel("Desplazamiento [in]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)");
    grid on;

    zitta = [0.1 0.1 0.1]; % relacion de amortiguamiento

    x2 = Descomposicion_Modal_Amortiguado(zitta, M, K, x0, dx0, t); % Calculo de la respuesta en el tiempo
    
    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo con amortiguamiento');
    hold on;
    plot(t,x2(1,:),"LineWidth",1.5,"Color","b");
    plot(t,x2(2,:),"LineWidth",1.5,"Color","r");
    plot(t,x2(3,:),"LineWidth",1.5,"Color","g");
    hold off;
    ylabel("Desplazamiento [in]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)");
    grid on;

    % Autovalores y autovectores
    [V,lambda] = eig(K,M);

    % Frecuencias naturales y amortiguadas
    wn = sqrt(diag(lambda)); % frecuencias naturales [rad/s] de cada modo de vibracion
    t_ev = 2*pi/wn(1) % periodo de evaluacion [s]
    X1 = [interp1(t,x1(1,:),t_ev),interp1(t,x1(2,:),t_ev),interp1(t,x1(3,:),t_ev)]; % desplazamientos para cada grado de libertad en el tiempo t_ev [s] para el caso sin amortiguamiento
    X2 = [interp1(t,x2(1,:),t_ev),interp1(t,x2(2,:),t_ev),interp1(t,x2(3,:),t_ev)]; % desplazamientos para cada grado de libertad en el tiempo t_ev [s] para el caso con amortiguamiento
    % Tabla que compara los desplazamientos para cada grado de libertad en el tiempo t_ev [s] para ambos casos
    T = table(X1',X2','VariableNames',{'Desplazamiento_No_Amortiguado','Desplazamiento_Amortiguado'},'RowNames',{'x1','x2','x3'});
    disp(T);


end

function x = Descomposicion_Modal_No_Amortiguado(M, K, x0, dx0, t)
    % Autovalores y autovectores
    [V,lambda] = eig(K,M);
    % disp("Autovalores:");disp(lambda);
    % disp("Autovectores:");disp(V);

    % disp('Las frecuencias angulares naturales en [rad/s] son:');
    wn = diag(sqrt(lambda)) % frecuencias naturales [rad/s] de cada modo de vibracion

    % Calculo de las frecuencias naturales en [Hz]
    % disp('Las frecuencias naturales en [Hz] son:');
    fn = wn/(2*pi); % frecuencias naturales [Hz] de cada modo de vibracion

    % Modos de vibracion normalizados al primer modo
    % disp('Modos de vibracion normalizados al primer modo:');
    V_norm = V./V(1,:) % modo de vibracion normalizado al primer modo

    % Relación de amortiguamiento nula
    zitta = 0; % coeficiente de amortiguamiento

    % Metodo de Descomposicion Modal

    M_modal = V'*M*V; % matriz de masa modal
    K_modal = V'*K*V; % matriz de rigidez modal

    % Condiciones iniciales en el espacio modal
    y0 = V' * M * x0; % desplazamiento inicial en el espacio modal
    dy0 = V' * M * dx0; % velocidad inicial en el espacio modal

    % Resolucion de las ecuaciones de movimiento en el espacio modal
    Y1 = zeros(length(wn),length(t)); % desplazamiento en el espacio modal

    for i = 1:length(wn)
        Y1(i,:) = y0(i)*cos(wn(i)*t) + dy0(i)/(wn(i))*sin(wn(i)*t);
    end

    % Paso de coordenadas del espacio modal al espacio fisico
    x = V*Y1; % desplazamiento en el espacio fisico
    
    
end

function x = Descomposicion_Modal_Amortiguado(zitta, M, K, x0, dx0, t)

    % Autovalores y autovectores
    [V,lambda] = eig(K,M);

    % Frecuencias naturales y amortiguadas
    wn = sqrt(diag(lambda)); % frecuencias naturales [rad/s] de cada modo de vibracion

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