function ej4_tp7;clc;close all;
    % parametros del sistema con 2 grados de libertad sin amortiguamiento
    m1 = 100; % masa [kg]
    m2 = m1/1000; % masa [kg]

    k1 = 1000; % rigidez [N/m]
    k2 = 3000; % rigidez [N/m]

    M = diag([m1,m2]); % matriz de masa
    K = [
          k1+k2  ,  -k2  ;
         -k2  , k2;
        ]; % matriz de rigidez

    zitta = [0 0]; % relacion de amortiguamiento

    % autovalores y autovectores
    [V,lambda] = eig(K,M);
    wn =diag(sqrt(lambda)) % frecuencias naturales [rad/s]
    % V_N = V./V(1,:) % autovectores normalizados

    % Condiciones iniciales
    x0 = [0 0]'; % desplazamiento inicial [m]
    dx0 = [0 0]'; % velocidad inicial [m/s]

    % Tiempo de analisis
    tf = 2; % tiempo final [s]
    dt = 1e-6; % paso de tiempo [s]
    t = 0:dt:tf; % vector de tiempo [s]

    % fuerza externa
    Po = 5000; % [N]
    wp = sqrt(k2/m2); % frecuencia de la carga [rad/s]

    P =zeros(length(M),length(t));
    P(1,:) = Po.*sin(wp*t); % fuerza externa en el tiempo

    [phi,rho,X] = Descomposicion_Modal_Amortiguado( zitta, M, K,V,wn, x0, dx0, t ,P,Po,wp);

    x1_max_teorico = 0;
    x2_max_teorico = Po/k2;
    x1_max = max(abs(X(1,:)));
    x2_max = max(abs(X(2,:)));

    % tabla para comparar los teoricos con los obtenidos
    T = table([x1_max_teorico; x2_max_teorico], [x1_max; x2_max], 'VariableNames', {'Teorico', 'Obtenido'}, 'RowNames', {'max(x1)', 'max(x2)'})

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo sin amortiguamiento');
    hold on;
    plot(t,X(1,:),"LineWidth",1.5,"Color","r");
    plot(t,X(2,:),"LineWidth",1.5,"Color","g");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)");
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
    PoModal = V' * [Po; 0];          %Paso las amplitudes de las fuerzas a coordenadas modales

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