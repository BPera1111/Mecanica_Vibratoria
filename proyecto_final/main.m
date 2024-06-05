function main;clc;close all; g = false; % graph flag

    if g; vel_power() ; end

    %% DEFINICIÓN DE PARÁMETROS DEL SISTEMA
    m1=26067; %% Masa 1
    m2=129209; %% Masa 2 
    
    z=0.015; zitta = [z z z z]; % relacion de amortiguamiento
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
    
    dt = 0.005; % time step
    t=0:dt:50; % time vector
    
    %% ENCUENTRO AUTOVECTORES Y AUTOVALORES
    [V,lambda]=eig(K,M); %el X que me larga ya esta normalizado de modo q X'*m*X=I
    
    wn=diag(sqrt(lambda)); %wn es un vector con las frecuencias naturales
    wd=wn.*sqrt(1-z.^2); %wd es un vector con las frecuencias de amortiguamiento
    
   

    %% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA
    FW=agua(t,g);
    
    %% DEFINO LA VELOCIDAD DEL VIENTO Y PARAMETROS DEL VIENTO
    % vel_viento=    ; %% Velocidad del viento [nudos]
    Fv = viento(t,g);
   
    %% DEFINO PARAMETROS Y FUERZA DE DESBALANCE

    Fd = desbalance(t,g);

    %% DEFINO LA FUERZA EXTERNA
    P = fuerza_externa(t,x0,FW,Fv,Fd);

    %% SOLUCIÓN DE LA ECUACIÓN DIFERENCIAL
    
    X = Descomposicion_Modal_Amortiguado(M,V,wn,t,P,wd,zitta,dt );
    trans = Descomposicion_Modal_Amortiguado_Transitoria( M,zitta,V,wn, x0, dx0, t);

    res = trans + X;

    res(:,40/dt)

    
    % respuesta(t,X,"Respuesta Permanente en el tiempo");

    % respuesta(t,trans,"Respuesta Transitoria en el tiempo");

    % respuesta(t,res,"Respuesta Total en el tiempo");
    

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

function FW=agua(t,graph) %% Función que devuelve la velocidad y la fuerza del agua en función del tiempo	
    %% DEFINO LA VELOCIDAD Y LA CARGA DEL AGUA

    Twa=5; %% Periodo del oleaje
    VAmax=2; %% Velocidad máxima del agua
    Cd=0.7; %% Coeficiente de arrastre
    RoW=1000; %% Densidad del agua
    Ainmersa=4*4.176; %% Área a considerar, h sobre agua 14m y diametro 4.176m
    VAt=zeros(1,length(t)); %% Vector agua
    for i=1:length(t)
       tmod=mod(t(i),Twa); %llevo el valor sub i de t a un valor dentro del periodo (entre 0 y Tp)
       VAt(i)= (tmod>=0).*(tmod<=Twa*4/5).*(((VAmax/2)/(Twa*4/5))*tmod+VAmax/2)+(tmod>Twa*4/5).*(tmod<=Twa).*(-VAmax/2*tmod+6);
    end

    VAt=VAt+VAmax/2;
    FW=0.5*RoW*Cd*Ainmersa.*VAt.^2.*(1-exp(-t./5));

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

function Fviento=viento(t,graph) %% Función que devuelve la fuerza del viento en función del tiempo
    
    % Agregar el calculo de la fuerza del viento en función del viento
    frecviento=1;
    Aviento=2000;
    for i=1:length(t)
       Fviento(i,1)=(100000+Aviento*sin(2*pi*frecviento*t(i))).*(1-exp(-t(i)/5));
    end

    if graph
        figure(1);
        subplot(4,1,3)
        plot(t,Fviento,'LineWidth',2)
        title('Fuerza del viento [N]')
    end
end

function Fdesb=desbalance(t,graph) %% Función que devuelve la fuerza de desbalance en función del tiempo

    F0Desbalance=1028.93;
    wRotor=2.572;

    Fdesb=F0Desbalance.*sin(wRotor*t);

    if graph
        figure(1);
        subplot(4,1,4)
        plot(t,Fdesb,'LineWidth',2)
        title('Fuerza de desbalance [N]')
    end

end

function P=fuerza_externa(t,x0,Fa,Fv,Fd) %% Función que devuelve la fuerza externa en función del tiempo


    P=zeros(length(x0),length(t)); % fuerza externa en el tiempo
    P(1,:) = 0; 
    P(2,:) = Fa;
    P(3,:) = Fd;
    P(4,:) = Fv;


end

function X = Descomposicion_Modal_Amortiguado(M,V,wn, t,P,wd,zitta,dt )
    % ----------------- MATRICES MODALES -----------------------------

    Mmodal = round(V' * M * V);
    Fmodal = V' * P;

    Y = zeros(length(wn), length(t));
    for i = 1:length(wn)
        for j = 1:length(t)
            integral_value = 0;
            for k = 1:j
                h = (1/Mmodal(i,i)) * exp(-zitta(i) * wn(i) * (t(j) - t(k))) * sin(wd(i) * (t(j) - t(k)));
                integral_value = integral_value + h * Fmodal(i,k) * dt;
            end
            Y(i,j) = integral_value;
        end
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
    plot(t,X(1,:),"b");
    plot(t,X(2,:),"r");
    plot(t,X(3,:),"g");
    plot(t,X(4,:),"k");
    hold off;
    ylabel("Desplazamiento [m]");
    xlabel("Tiempo [s]");
    legend("x1(t)","x2(t)","x3(t)","x4(t)");
    grid on;
    
end

