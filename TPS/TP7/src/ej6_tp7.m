function ej6_tp7;clc;close all;
    %Masas (seleccionadas de forma arbitraria, pero teniendo en cuenta que)
    m=3000;
    m1=m;
    m2=4*m;
    m3=m;
    M = diag([m1,m2,m3]);

    %Rigidez
    E=6.9*10^9;
    I=5.2*10^(-6);
    l=2;
    k=3*E*I/l;

    K = [k    -k    0;
        -k    2*k  -k;
        0    -k    k];

    %Condiciones iniciales
    x0 = [0.2; 0; 0]; %Desplazamientos iniciales
    dx0 = [0; 0; 0] ; %Velocidades iniciales

    %Tiempo para graficar
    tf=10;
    dt=0.01;
    t=0:dt:tf;

    %Relacion de amortiguamiento crítico
    zitta = [0; 0; 0];

    %% ------------------ AUTOVECTORES Y AUTOVALORES ------------------

    [V,lambda] = eig(K,M);
    w = diag(sqrt(lambda))
    V_norm = V ./ V(1,:)

    %% ------------------ RESPUESTA EN EL TIEMPO ------------------

    x = Descomposicion_Modal_Amortiguado(zitta, M,V,w, x0, dx0, t);

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo con amortiguamiento');
    hold on;
    plot(t,x(1,:),"LineWidth",1.5,"Color","b");
    plot(t,x(2,:),"LineWidth",1.5,"Color","r");
    plot(t,x(3,:),"LineWidth",1.5,"Color","g");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)");
    grid on;

    % Tabala para mostrar los maximos desplazamientos y sus tiempos
    maximos = [max(x(1,:)), max(x(2,:)), max(x(3,:))];
    tiempos = [t(find(x(1,:)==maximos(1))), t(find(x(2,:)==maximos(2))), t(find(x(3,:)==maximos(3)))];
    T = table(maximos', tiempos', 'VariableNames', {'DesplazamientoMaximo', 'TiempoDesplazamientoMaximo'},'RowNames', {'max(x1)', 'max(x2)','max(x3)'})

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