function ej2_tp7;clc;close all;
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

    zitta = [0.1 0.1 0.1]; % relacion de amortiguamiento

    % autovalores y autovectores
    [V,lambda] = eig(K,M);
    wn =diag(sqrt(lambda)) % frecuencias naturales [rad/s]
    V_N = V./V(1,:) % autovectores normalizados

    

    % Condiciones iniciales
    x0 = [ 0.3; -0.8 ; 0.3 ]; % desplazamiento inicial [in]
    dx0 = [ 0; 0; 0 ]; % velocidad inicial [in/s]

    % Tiempo de analisis
    tf = 2; % tiempo final [s]
    dt = 1e-6; % paso de tiempo [s]
    t = 0:dt:tf; % vector de tiempo [s]

    % fuerza externa
    Po_N = 5000; % [N]
    Po = Po_N/4.44822; % [kips]
    wp = 1.1*wn(1); % frecuencia de la carga [rad/s]

    P =zeros(length(M),length(t));
    P(1,:) = Po.*sin(wp*t); % fuerza externa en el tiempo

    [phi,rho,X] = Descomposicion_Modal_Amortiguado( zitta, M, K,V,wn, x0, dx0, t ,P,Po,wp);

    phi
    rho

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo con amortiguamiento');
    hold on;
    plot(t,X(1,:),"LineWidth",1.5,"Color","b");
    plot(t,X(2,:),"LineWidth",1.5,"Color","r");
    plot(t,X(3,:),"LineWidth",1.5,"Color","g");
    hold off;
    ylabel("Desplazamiento [in]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)");
    grid on;


end

function [phi,rho,X] = Descomposicion_Modal_Amortiguado( zitta, M, K,V,w, x0, dx0, t,P,Po,wp )
    % ----------------- MATRICES MODALES -----------------------------

    Mmodal = round(V' * M * V);
    Kmodal = round(V' * K * V);
    Cmodal = 2*zitta.*w.*eye(size(M));
    Y0 = V' * M * x0;                   %desplazamiento inicial en coordenadas modales
    dY0 = V' * M * dx0;                 %velocidad inicial en coordenadas modales
    Pmodal = V' * P;                    %Paso fuerza a coordenadas modales
    PoModal = V' * [Po; 0; 0];          %Paso las amplitudes de las fuerzas a coordenadas modales

    % ----------------- SOLUCION A LA EDO ----------------------------

    e=exp(1);
    wd = w.*sqrt(1-zitta.^2); %frecuencia natural del sistema amortiguado
    beta = zeros(length(w),1);
    phi = zeros(length(w),1);
    rho = zeros(length(w),1);
    H = zeros(length(w),1);


    for i=1:length(w)
        beta(i)=wp/w(i);
        H(i) = 1/sqrt( (1-beta(i)^2)^2 + (2*zitta(i)*beta(i))^2 );
        phi(i)=atan(2*zitta(i)*beta(i)/(1-beta(i)^2));
        PoModal(i)=max(abs(Pmodal(i,:)));
        rho(i)=(PoModal(i)/Kmodal(i,i))*H(i);
    end

    Y = zeros(length(w), length(t));

    for i = 1:length(w)    
        Y(i,:) = rho(i,:).*sin(wp*t-phi(i));
    end
    %Ahora paso de coordenadas modales a geométricas
    X = V * Y;
    
end