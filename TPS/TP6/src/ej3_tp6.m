function ej3_tp6;clc;close all;
    % Parametros del sistema
    m = 0.1; % Masa [kips.sec^2/in]
    k = 5; % Rigidez [kips/in]
    c = 0.2; % Amortiguamiento [kips.sec/in]
    wn = sqrt(k/m); % Frecuencia natural [rad/s]
    % zitta = c/(2*sqrt(m*k)); % Coeficiente de amortiguamiento
    % wd = wn*sqrt(1-zitta^2); % Frecuencia de amortiguamiento [rad/s]

    dt = 0.1; % Paso de tiempo [s]
    tf = 8; % Tiempo final [s]
    t = 0:dt:tf; % Vector de tiempo [s]

    % Parametros de la carga externa
    P = carga(t); % Carga externa en el tiempo

    [x,dx] = Newmark_Beta(wn,t,dt,P,m,c,k); % Calculo la respuesta en el tiempo con el metodo de Newmark Beta

    % Fuerza del resorte
    Fk = (-k)*x;

    % Mostramos con subplot la carga, el desplazamiento y la fuerza del resorte
    figure('Name', 'Respuesta en el tiempo');
    subplot(3,1,1);
    hold on;
    plot(t,P,"LineWidth",1.5,"Color","r");
    scatter(t,P,"filled",'r');
    hold off;
    ylabel("Carga [kips]");
    grid on;
    grid minor;

    subplot(3,1,2);
    hold on;
    plot(t,x,"LineWidth",1.5,"Color","b");
    scatter(t,x,"filled",'b');
    hold off;
    ylabel("Desplazamiento [in]");
    grid on;
    grid minor;

    subplot(3,1,3);
    hold on;
    plot(t,Fk,"LineWidth",1.5,"Color","g");
    scatter(t,Fk,"filled",'g');
    hold off;
    ylabel("Fuerza Elastica [kips]");
    xlabel("Tiempo [s]");
    grid on;
    grid minor;

    % mostar una tabla con los primeros 10 valores de t, P, x y Fk
    n = 0:1:10;
    T = table(n(1:10)',t(1:10)',P(1:10)',x(1:10)',dx(1:10)',Fk(1:10)','VariableNames',{'n','t','Pn','x(t)','v(t)','Fk(t)'});
    disp(T);
    


end

function P=carga(t) % Carga externa en el tiempo
    t_app = 0:0.1:0.8;
    P_app = [0,5,8,7,5,3,2,1,0];
    P = interp1(t_app,P_app,t,"linear",0);
end

function [x,dx,ddx] = Newmark_Beta(T,t,dt,P,m,c,k) % Metodo de Newmark Beta para coef de amortiguamiento constante
    if dt/T<sqrt(3)/pi
        %Inicializo vectores de desplazamiento, velocidad y aceleracion
        x=zeros(1,length(t));
        dx=zeros(1,length(t));
        ddx=zeros(1,length(t));
        %Establezco condiciones iniciales
        x(1)=0;
        dx(1)=0;
        ddx(1)=0;

        %Calculo la rigidez efectiva
        keff=k+3*c/(dt)+6*m/(dt^2);

        %inicializo el vector de carga efectiva
        Peff=zeros(1,length(t));

        %Ecuaciones de Newmark Beta
        for i=2:length(t)
            Peff(i)=P(i)+m*(6*x(i-1)/dt^2 + 6*dx(i-1)/dt + 2*ddx(i-1))+c*(3*x(i-1)/dt + 2*dx(i-1)+dt*ddx(i-1)/2);
            x(i)=1/keff*Peff(i);
            dx(i) = -2*dx(i-1) -dt*ddx(i-1)/2 + 3/dt*(x(i)-x(i-1));
            ddx(i)=6/dt^2*(x(i)-x(i-1))-6/dt*dx(i-1)-2*ddx(i-1);
        end
    else
        disp("Error: El paso de tiempo es mayor al permitido por el metodo de Newmark Beta");
    end
    
end