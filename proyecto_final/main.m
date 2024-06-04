function main;clc;close all; g = false; % graph flag

    if g; vel_power() ; end

    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    m1=26067; %% Masa 1
    m2=129209; %% Masa 2 
    
    z=0.015;
    zitta = [z z z z]; % relacion de amortiguamiento
    x0 = [ 0 ; 0 ; 0 ; 0 ]; % desplazamiento inicial [m]
    dx0 = [ 0; 0 ; 0 ; 0 ]; % velocidad inicial [m/s]
    
    k1=2.01e8; %% Rigidez 1 [N/m]
    k2=9.36e6; %% Rigidez 2 [N/m]
    
    M=[m1 0 0 0;
        0 m1 0 0;
        0 0 m2 0;
        0 0 0 m2]; %% Matriz de masa
    
    K=[k1+k2 0 -k2 0;
        0 k1+k2 0 -k2;
        -k2 0   k2 0;
         0 -k2  0 k2]; %% Matriz de Rigidez
    
    dt = 0.001; % time step
    t=0:dt:10; % time vector
    
    %% ENCUENTRO AUTOVECTORES Y AUTOVALORES
    [V,lambda]=eig(K,M); %el X que me larga ya esta normalizado de modo q X'*m*X=I
    
    wn=diag(sqrt(lambda)); %wn es un vector con las frecuencias naturales
    wd=wn.*sqrt(1-z.^2); %wd es un vector con las frecuencias de amortiguamiento
    
   
    %% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA
    
    Twa=5; %% Periodo del oleaje
    VAmax=3; %% Velocidad máxima del agua
    Cd=0.7; %% Coeficiente de arrastre
    RoW=1000; %% Densidad del agua
    Ainmersa=14*4.176; %% Área a considerar, h sobre agua 14m y diametro 4.176m
    
    [VAt,FW]=agua(t,Twa,VAmax,Cd,RoW,Ainmersa,g);
    
    %% DEFINO LA VELOCIDAD DEL VIENTO Y PARAMETROS DEL VIENTO
    % vel_viento=    ; %% Velocidad del viento [nudos]
    Fv = viento(t,g);
   
    %% DEFINO PARAMETROS Y FUERZA DE DESBALANCE

    F0Desbalance=1028.93;
    wRotor=2.572;

    Fd = desbalance(t,F0Desbalance,wRotor,g);

    %% DEFINO LA FUERZA EXTERNA
    P = fuerza_externa(t,x0,FW,Fv,Fd);


    %% SOLUCIÓN DE LA ECUACIÓN DIFERENCIAL
    wp = [ 2*pi/Twa 0 0 wRotor] % frecuencia de la carga [rad/s]
    Po = [ max(FW) 0 Fv(1) F0Desbalance]     % Amplitudes de las fuerzas [N]
    
    [phi,rho,X] = Descomposicion_Modal_Amortiguado( zitta, M, K,V,wn, x0, dx0, t,P,Po,wp );
    trans = Descomposicion_Modal_Amortiguado_Transitoria( M,zitta,V,wn, x0, dx0, t);

    res = trans + X;

    phi;
    rho;
    
    respuesta(t,X,"Respuesta Permanente en el tiempo");

    respuesta(t,trans,"Respuesta Transitoria en el tiempo");

    respuesta(t,res,"Respuesta Total en el tiempo");
    

end

function vel_power % Power vs Velocity
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

function [VAt,FW]=agua(t,Twa,VAmax,Cd,RoW,Ainmersa,graph) %% Función que devuelve la velocidad y la fuerza del agua en función del tiempo	
    VAt=zeros(1,length(t)); %% Vector agua
    for i=1:length(t)

        tmod=mod(t(i),Twa); %llevo el valor sub i de t a un valor dentro del periodo (entre 0 y Tp)
        VAt(i)= (tmod>=0).*(tmod<=Twa).*(((VAmax/2)/Twa)*tmod);
        
    end

    VAt=VAt+VAmax/2;
    FW=0.5*RoW*Cd*Ainmersa.*VAt.^2;

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

function [Fviento]=viento(t,graph) %% Función que devuelve la fuerza del viento en función del tiempo
    
    % Agregar el calculo de la fuerza del viento en función del viento

    for i=1:length(t)
        Fviento(i)=101213.24; %% Fviento 33kN en todo T
    end

    if graph
        figure(1);
        subplot(4,1,3)
        plot(t,Fviento,'LineWidth',2)
        title('Fuerza del viento [N]')
    end
end

function [Fdesb]=desbalance(t,F0Desbalance,wRotor,graph) %% Función que devuelve la fuerza de desbalance en función del tiempo

    Fdesb=F0Desbalance.*sin(wRotor*t);

    if graph
        figure(1);
        subplot(4,1,4)
        plot(t,Fdesb,'LineWidth',2)
        title('Fuerza de desbalance [N]')
    end

end

function [P]=fuerza_externa(t,x0,Fa,Fv,Fd) %% Función que devuelve la fuerza externa en función del tiempo

    P=zeros(length(x0),length(t)); % fuerza externa en el tiempo
    P(1,:) = 0; 
    P(2,:) = Fa;
    P(3,:) = Fd;
    P(4,:) = Fv;


end

function [phi,rho,X] = Descomposicion_Modal_Amortiguado( zitta, M, K,V,w, x0, dx0, t,P,Po,wp )
    % ----------------- MATRICES MODALES -----------------------------

    % Mmodal = round(V' * M * V);
    Kmodal = round(V' * K * V);
    % Cmodal = 2*zitta.*w.*eye(size(M));
    Y0 = V' * M * x0;                   %desplazamiento inicial en coordenadas modales
    dY0 = V' * M * dx0;                 %velocidad inicial en coordenadas modales
    Pmodal = V' * P;                    %Paso fuerza a coordenadas modales
    PoModal = V' .* Po;          %Paso las amplitudes de las fuerzas a coordenadas modales

    % ----------------- SOLUCION A LA EDO ----------------------------

    % e=exp(1);
    % wd = w.*sqrt(1-zitta.^2); %frecuencia natural del sistema amortiguado
    beta = zeros(length(w),1);
    phi = zeros(length(w),1);
    rho = zeros(length(w),1);
    H = zeros(length(w),1);


    for i=1:length(w)
        beta(i)=wp(i)/w(i);
        H(i) = 1/sqrt( (1-beta(i)^2)^2 + (2*zitta(i)*beta(i))^2 );
        phi(i)=atan(2*zitta(i)*beta(i)/(1-beta(i)^2));
        PoModal(i)=max(abs(Pmodal(i,:)));
        rho(i)=(PoModal(i)/Kmodal(i,i))*H(i);
    end

    Y = zeros(length(w), length(t));

    for i = 1:length(w)    
        Y(i,:) = rho(i,:).*sin(wp(i)*t-phi(i));
    end
    %Ahora paso de coordenadas modales a geométricas
    X = V * Y ;
    
end

function x = Descomposicion_Modal_Amortiguado_Transitoria( M,zitta,V,wn, x0, dx0, t)

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

function respuesta(t,X,titulo)
    % Mostramos la respuesta en el tiempo
    figure('Name', titulo);
    hold on;
    plot(t,X(1,:),"LineWidth",1.5,"Color","b");
    plot(t,X(2,:),"LineWidth",1.5,"Color","r");
    plot(t,X(3,:),"LineWidth",1.5,"Color","g");
    plot(t,X(4,:),"LineWidth",1.5,"Color","y");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)","x4(t)");
    grid on;
    
end