function main;clc;close all; g = true; % graph flag

    if g; vel_power() ; end %#ok<UNRCH>

    
    dt = 0.005; % time step
    t=0:dt:50; % time vector
    
    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    m1=62000; %% Masa 1
    m2=93000; %% Masa 2
    m3=0.03*(m1+m2); %% Masa TMD

    z=[0.015;0.015;0.05];

    k1=2e8; %% Rigidez 1
    k2=9e6; %% Rigidez 2
    k3=(9.4)^2*m3; %% Rigidez 3

    M= [m1 0 0;
        0 m2 0;
        0 0 m3] %% Matriz de masa

    K=[k1+k2 -k2 0;
        -k2   k2+k3 -k3;
        0 -k3 k3] %% Matriz de Rigidez

    xINI=[0;0;0];
    dxINI=[0;0;0];

    %% ENCUENTRO AUTOVECTORES Y AUTOVALORES
    [V,lambda]=eig(K,M); %el X que me larga ya esta normalizado de modo q X'*m*X=I

    wn=diag(sqrt(lambda));
    wd=wn.*sqrt(1-z.^2);
    

    %% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA
    FW=agua(t,g);
    
    %% DEFINO LA VELOCIDAD DEL VIENTO Y PARAMETROS DEL VIENTO
    vel_viento=7.64;    ; %% Velocidad del viento [m/s]
    Fv = viento(t,vel_viento,g);
   
    %% DEFINO PARAMETROS Y FUERZA DE DESBALANCE

    Fd = desbalance(t,vel_viento,g);

    %% EJECUTO EL SISTEMA

    ejecutar(FW,Fv,Fd,t,dt);

end

function vel_power %#ok<DEFNU> % Power vs Velocity
    vel = 0:0.1:27; % velocity vector
    power = zeros(1, length(vel)); % power vector
    v_cin = 3;
    v_rat = 12;
    v_cou = 24;
    p_r = 3;

    q =  p_r*((vel.^2-v_cin^2)/(v_rat^2-v_cin^2));
    
    % set power values based on velocity
    power(vel <= v_cin) = 0;
    power(vel > v_cin & vel <= v_rat) = q(vel > v_cin & vel <= v_rat);
    power(vel > v_rat & vel <= v_cou) = p_r;
    power(vel > v_cou) = 0;
    
    % plot power vs velocity
    figure(2);
    plot(vel, power, 'LineWidth', 2);
    xlabel('Velocity [m/s]');
    ylabel('Power[MW]');
    % legend('Power vs Velocity');
    title('Power vs Velocity');
end

function FW=agua(t,graph) %% Función que devuelve la velocidad y la fuerza del agua en función del tiempo	
    %% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA

    Twa = 5; %% Periodo del oleaje
    VAmax = 2; %% Velocidad máxima del agua
    Cd = 0.7; %% Coeficiente de arrastre
    RoW = 1000; %% Densidad del agua
    Ainmersa = 14 * 4.176; %% Área a considerar, h sobre agua 14m y diametro 4.176m

    VAt = zeros(1, length(t)); %% Vector agua
    for i = 1:length(t)

        tmod = mod(t(i), Twa); % llevo el valor sub i de t a un valor dentro del periodo (entre 0 y Tp)
        VAt(i) = (tmod >= 0) .* (tmod <= Twa * 4 / 5) .* (((VAmax / 2) / (Twa * 4 / 5)) * tmod + VAmax / 2) + (tmod > Twa * 4 / 5) .* (tmod <= Twa) .* (-VAmax / 2 * tmod + 6);
    end


    FW = (0.5 * RoW * Cd * Ainmersa .* VAt.^2 + 91120);
    

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

function Fviento=viento(t,V_viento,graph) %% Función que devuelve la fuerza del viento en función del tiempo
    
    Cd=0.7;
    P_aire=1.203;
    D=4.176;
    A_mastil=(84-14)*D;
    
    Fv_mastil=(1/2)*Cd*P_aire*A_mastil*V_viento^2;
    
    Fv_aspas=(1/2)*Cd*P_aire*(0.4*pi*56^2)*V_viento^2;
    
    Avientob=Fv_mastil + Fv_aspas;


    Aviento=20000;
    % Avientob=100000;
    frecViento=0.5;

    %Fviento = 100000 * ones(1, length(t));

    %Agrego variacion al viento
    for i=1:2500
        Fviento(i)=Avientob+Aviento*sin(2*pi*frecViento*t(i))+Aviento/3*sin(3*2*pi*frecViento*t(i))+Aviento/4*cos(4.1*2*pi*frecViento*t(i))+Aviento/8.2*sin(8.34*2*pi*frecViento*t(i)); %#ok<AGROW>
    end
    for i=2501:5000
        Fviento(i)=Avientob+Aviento/1.2*sin(1.1*2*pi*frecViento*t(i))+Aviento/3.23*sin(4.21*2*pi*frecViento*t(i))+Aviento/5.23*cos(6.23*2*pi*frecViento*t(i))+Aviento/8.26*sin(5.32*2*pi*frecViento*t(i));
    end
    for i=5001:7500
        Fviento(i)=Avientob+Aviento/1.77*cos(1.1*2*pi*frecViento*t(i))+Aviento/3.91*sin(2.21*2*pi*frecViento*t(i))+Aviento/3.678*cos(7.1*2*pi*frecViento*t(i))+Aviento/9.4*sin(4.375*2*pi*frecViento*t(i));
    end
    for i=7501:10001
        Fviento(i)=Avientob+Aviento/0.992*cos(1.632*2*pi*frecViento*t(i))+Aviento/5.13*sin(4.581*2*pi*frecViento*t(i))+Aviento/6.8*sin(7.1*2*pi*frecViento*t(i))+Aviento/9.4*sin(4.65*2*pi*frecViento*t(i));
    end

    if graph
        figure(1);
        subplot(4,1,3)
        plot(t,Fviento,'LineWidth',2)
        title('Fuerza del viento [N]')
    end
end

function Fdesb=desbalance(t,V_viento,graph) %% Función que devuelve la fuerza de desbalance en función del tiempo

    R = 56;
    Uper=155.7;
    TSR=7;
    wRotor=(TSR/R)*V_viento;
    F0Desbalance=Uper*wRotor^2;
    Fdesb=F0Desbalance.*sin(wRotor*t);

    if graph
        figure(1);
        subplot(4,1,4)
        plot(t,Fdesb,'LineWidth',2)
        title('Fuerza de desbalance [N]')
    end

end

function Py=Fuerza_externa_y(t,x0,Fd) %% Función que devuelve la fuerza externa en función del tiempo
    Py=zeros(length(x0),length(t)); % fuerza externa en el tiempo
    Py(2,:) = Fd; 
end

function Px=Fuerza_externa_x(t,x0,Fa,Fv) %% Función que devuelve la fuerza externa en función del tiempo
    Px=zeros(length(x0),length(t)); % fuerza externa en el tiempo
    Px(1,:) = Fa;
    Px(2,:) = Fv;
end

function X = Fact_Magnificacion_Dinam(M,K,wn,wd,P,X,t)
    % Transformar fuerzas a coordenadas modales
    Mmodal=round(X'*M*X); %M modal
    Kmodal=round(X'*K*X); %K modal
    Fmodal = X' * P;

    % Respuesta permanente modal
    Y_perm = zeros(length(wn), length(t));
    for i = 1:length(wn)
        for j = 1:length(t)
            Y_perm(i,j) = (1/(Kmodal(i,i) - Mmodal(i,i) * wd(i)^2)) * Fmodal(i,j);
        end
    end

    % Paso de coordenadas modales a geométricas
    X = X * Y_perm;

end

function X = Transformada_De_Fourier(M,K,t,P,dt,V,wn,z)
    Mmodal = round(V' * M * V); % M modal
    Kmodal = round(V' * K * V); % K modal
    Cmodal = 2*z.*wn.*eye(size(Mmodal)); % C modal

    % Transformar fuerzas a coordenadas modales
    Fmodal = V' * P;

    % Frecuencias de Fourier
    omega = 2 * pi * (0:(length(t)-1)) / (length(t) * dt);

    % Función de transferencia H(iw)
    H = zeros(size(Mmodal, 1), length(omega));
    for i = 1:length(omega)
        H(:, i) = diag(((-omega(i)^2 * Mmodal + 1i * omega(i) * Cmodal + Kmodal) \ eye(size(Mmodal))));
    end

    % Respuesta en frecuencia
    Y_freq = H .* fft(Fmodal, [], 2);

    % Transformar de vuelta al dominio del tiempo
    Y_perm = ifft(Y_freq, [], 2, 'symmetric');

    % Transformar respuesta permanente a coordenadas físicas
    X = V * Y_perm;

end

function x = Respuesta_Transitoria( M,zitta,V,wn, x0, dx0, t)

    % Paso de coordenadas del espacio fisico al espacio modal
    y0 = V'*M*x0; % desplazamiento inicial en el espacio modal
    dy0 = V'*M*dx0; % velocidad inicial en el espacio modal
    
    % C_modal = V'*C*V; % matriz de amortiguamiento modal [kips.sec/in]

    % zitta = diag(C_modal)./(2*wn); % relacion de amortiguamiento

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

function respuesta_1eje_tmd(t,X,titulo)
    % Mostramos la respuesta en el tiempo
    figure('Name', titulo);
    hold on;
    plot(t,X(1,:),"r");
    plot(t,X(2,:),"g");
    plot(t,X(3,:),"b");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","TMD(t)");
    grid on;
    
end

function respuesta_1eje(t,X,titulo)
    % Mostramos la respuesta en el tiempo
    figure('Name', titulo);
    hold on;
    plot(t,X(1,:),"r");
    plot(t,X(2,:),"g");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)");
    grid on;
    
end

function [M,K,wn,wd,z,xINI,dxINI,V]=eje_ind_tmd_forzadas()
    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    m1=62000; %% Masa 1
    m2=93000; %% Masa 2
    m3=0.03*(m1+m2); %% Masa TMD

    z=[0.015;0.015;0.05];

    k1=2e8; %% Rigidez 1
    k2=9e6; %% Rigidez 2
    k3=(9.4)^2*m3; %% Rigidez 3

    M= [m1 0 0;
        0 m2 0;
        0 0 m3]; %% Matriz de masa

    K=[k1+k2 -k2 0;
        -k2   k2+k3 -k3;
        0 -k3 k3]; %% Matriz de Rigidez

    xINI=[0;0;0];
    dxINI=[0;0;0];

    %% ENCUENTRO AUTOVECTORES Y AUTOVALORES
    [V,lambda]=eig(K,M); %el X que me larga ya esta normalizado de modo q X'*m*X=I

    wn=diag(sqrt(lambda));
    wd=wn.*sqrt(1-z.^2);

end

function [M,K,wn,wd,z,xINI,dxINI,V]=eje_ind_forzadas()
    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    m1=62000; %% Masa 1
    m2=93000; %% Masa 2
    %m3=0.05*(m1+m2); %% Masa TMD

    z=[0.015;0.015];

    k1=2e8; %% Rigidez 1
    k2=9e6; %% Rigidez 2
    %k3=(20)^2*m3; %% Rigidez 3

    M= [m1 0;
        0 m2]; %% Matriz de masa

    K=[k1+k2 -k2;
        -k2   k2]; %% Matriz de Rigidez


    xINI=[0;0];
    dxINI=[0;0];

    %% ENCUENTRO AUTOVECTORES Y AUTOVALORES
    [V,lambda]=eig(K,M); %el X que me larga ya esta normalizado de modo q X'*m*X=I

    wn=diag(sqrt(lambda));
    wd=wn.*sqrt(1-z.^2);

end

function A = amplitud(x,t1,t2,dt)
    A = max((x(t1/dt:t2/dt)));
end


function m=num_random(m)
    m(2,2)=9e6;
    x1= m(2,2)*rand(); %generamos un número aleatorio entre 0 - 1 multiplicado por el valor para poder generar la carga aleatoria del viento
end

function ejecutar(FW,Fv,Fd,t,dt)

    %% EJE Y CON TMD FORZADO
    [M_yc, K_yc, wn_yc, wd_yc, zitta_yc, x0_yc, dx0_yc, V_yc] = eje_ind_tmd_forzadas();
    K_yc = num_random(K_yc);
    P_Yc = Fuerza_externa_y(t,x0_yc,Fd);

    y_transitoria_c = Respuesta_Transitoria( M_yc,zitta_yc,V_yc,wn_yc, x0_yc, dx0_yc, t);
    y_permanente_2c = Fact_Magnificacion_Dinam(M_yc,K_yc,wn_yc,wd_yc,P_Yc,V_yc,t);

    y_res_2c = y_transitoria_c + y_permanente_2c;

    respuesta_1eje_tmd(t,y_res_2c,"EJE Y CON TMD 2");



    %% EJE Y SIN TMD FORZADO
    [M_ys, K_ys, wn_ys, wd_ys, zitta_ys, x0_ys, dx0_ys, V_ys] = eje_ind_forzadas();
    K_ys = num_random(K_ys);
    P_Ys = Fuerza_externa_y(t,x0_ys,Fd);

    y_transitoria_s = Respuesta_Transitoria( M_ys,zitta_ys,V_ys,wn_ys, x0_ys, dx0_ys, t);
    y_permanente_2s = Fact_Magnificacion_Dinam(M_ys,K_ys,wn_ys,wd_ys,P_Ys,V_ys,t);

    y_res_2s = y_transitoria_s + y_permanente_2s;


    respuesta_1eje(t,y_res_2s,"EJE Y SIN TMD 2");


    %% EJE X CON TMD FORZADO
    [M_xc, K_xc, wn_xc, wd_xc, zitta_xc, x0_xc, dx0_xc, V_xc] = eje_ind_tmd_forzadas();
    K_xc = num_random(K_xc);
    P_Xc = Fuerza_externa_x(t,x0_xc,FW,Fv);

    x_transitoria_c = Respuesta_Transitoria( M_xc,zitta_xc,V_xc,wn_xc, x0_xc, dx0_xc, t);
    x_permanente_2c = Transformada_De_Fourier(M_xc,K_xc,t,P_Xc,dt,V_xc,wn_xc,zitta_xc);

    x_res_2c = x_transitoria_c + x_permanente_2c;
    
    respuesta_1eje_tmd(t,x_res_2c,"EJE X CON TMD 2");


    %% EJE X SIN TMD FORZADO
    [M_xs, K_xs, wn_xs, wd_xs, zitta_xs, x0_xs, dx0_xs, V_xs] = eje_ind_forzadas();
    K_xs = num_random(K_xs);
    P_Xs = Fuerza_externa_x(t,x0_xs,FW,Fv);

    x_transitoria_s = Respuesta_Transitoria( M_xs,zitta_xs,V_xs,wn_xs, x0_xs, dx0_xs, t);
    x_permanente_2s = Transformada_De_Fourier(M_xs,K_xs,t,P_Xs,dt,V_xs,wn_xs,zitta_xs);

    x_res_2s = x_transitoria_s + x_permanente_2s;

    respuesta_1eje(t,x_res_2s,"EJE X SIN TMD 2");


    t1 = 2
    t2 = 18
    fprintf('Tabla de resultados:\n');
    fprintf('----------------------------------------\n');
    Amp_TMDy = amplitud(y_res_2c(2,:),t1,t2,dt);
    Amp_SIN_TMDy = amplitud(y_res_2s(2,:),t1,t2,dt);
    fprintf('Eje Y con TMD forzada 2: %f\n', Amp_TMDy);
    fprintf('Eje Y sin TMD forzada 2: %f\n', Amp_SIN_TMDy);    
    reduction_percentage = (Amp_TMDy / Amp_SIN_TMDy) * 100;
    fprintf('Porcentaje de reducción al aplicar el TMD: %.2f%%\n', reduction_percentage);
    fprintf('----------------------------------------\n');

    Amp_TMDx = amplitud((x_res_2c(2,:)),t1,t2,dt);
    Amp_SIN_TMDx = amplitud(x_res_2s(2,:),t1,t2,dt);
    fprintf('Eje X con TMD forzada 2: %f\n', Amp_TMDx);
    fprintf('Eje X sin TMD forzada 2: %f\n', Amp_SIN_TMDx);
    reduction_percentage = (Amp_TMDx / Amp_SIN_TMDx) * 100;
    fprintf('Porcentaje de reducción al aplicar el TMD: %.2f%%\n', reduction_percentage);
    fprintf('----------------------------------------\n');

    t1 = 30
    t2 = 50
    fprintf('Tabla de resultados:\n');
    fprintf('----------------------------------------\n');
    Amp_TMDy = amplitud(y_res_2c(2,:),t1,t2,dt);
    Amp_SIN_TMDy = amplitud(y_res_2s(2,:),t1,t2,dt);
    fprintf('Eje Y con TMD forzada 2: %f\n', Amp_TMDy);
    fprintf('Eje Y sin TMD forzada 2: %f\n', Amp_SIN_TMDy);    
    reduction_percentage = (Amp_TMDy / Amp_SIN_TMDy) * 100;
    fprintf('Porcentaje de reducción al aplicar el TMD: %.2f%%\n', reduction_percentage);
    fprintf('----------------------------------------\n');

    Amp_TMDx = amplitud((x_res_2c(2,:)),t1,t2,dt);
    Amp_SIN_TMDx = amplitud(x_res_2s(2,:),t1,t2,dt);
    fprintf('Eje X con TMD forzada 2: %f\n', Amp_TMDx);
    fprintf('Eje X sin TMD forzada 2: %f\n', Amp_SIN_TMDx);
    reduction_percentage = (Amp_TMDx / Amp_SIN_TMDx) * 100;
    fprintf('Porcentaje de reducción al aplicar el TMD: %.2f%%\n', reduction_percentage);
    fprintf('----------------------------------------\n');
    % % fprintf('Eje X e Y con TMD forzada 1: %f\n', res_1c(4,40/dt));
    % fprintf('Eje X e Y con TMD forzada 2: %f\n', res_2c(4,40/dt));
    % % % fprintf('Eje X e Y sin TMD forzada 1: %f\n', res_1s(4,40/dt));
    % fprintf('Eje X e Y sin TMD forzada 2: %f\n', res_2s(4,40/dt));
    % reduction_percentage_x4 = (res_2c(4,40/dt) / res_2s(4,40/dt)) * 100;
    % fprintf('Porcentaje de reducción al aplicar el TMD Viento Agua: %.2f%%\n', reduction_percentage_x4);
    % reduction_percentage_x3 = (res_2c(3,40/dt) / res_2s(3,40/dt)) * 100;
    % fprintf('Porcentaje de reducción al aplicar el TMD Desbalance: %.2f%%\n', reduction_percentage_x3);
    % fprintf('----------------------------------------\n');

    figure('Name', 'EJE x con TMD vs EJE x sin TMD');
    hold on;
    % plot(t,(x_res_2c(2,:)*0.6+(-mean(x_res_2c(2,:)*0.6)+mean(x_res_2s(2,:)))),"r");
    plot(t,x_res_2c(2,:),"r");
    plot(t,x_res_2s(2,:),"g");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("CON TMD","SIN TMD");
    grid on;
    figure('Name', 'EJE Y con TMD vs EJE Y sin TMD');
    hold on;
    plot(t,y_res_2c(2,:),"r");
    plot(t,y_res_2s(2,:),"g");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("CON TMD","SIN TMD");
    grid on;

end