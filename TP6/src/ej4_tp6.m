function ej4_tp6;clc;close all;
    % Parametros del sistema con un grado de libertad
    m = 0.2; % Masa [kips.sec^2/in]
    k = 8; % Rigidez [kips/in]
    c = 0.4; % Amortiguamiento [kips.sec/in]
    wn = sqrt(k/m); % Frecuencia natural [rad/s]
    % zitta = c/(2*sqrt(m*k)); % Coeficiente de amortiguamiento
    % wd = wn*sqrt(1-zitta^2); % Frecuencia de amortiguamiento [rad/s]
    T = 2*pi/wn; % Periodo [s]

    dt = 0.006; % Paso de tiempo [s]
    tf = 6; % Tiempo final [s]
    t = 0:dt:tf; % Vector de tiempo [s]

    % Parametros de la carga externa
    P = carga(t); % Carga externa en el tiempo

    [xt,dx] = PaP_dif_central(T,t,dt,P,m,c,k); % Calculo de la integral de Duhamel

    % fuerza del resorte
    Fk = (-k)*xt;

  

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo');
    yyaxis left;
    hold on;
    plot(t,xt,"LineWidth",1.5)
    % scatter(t,xt,"filled");  
    hold off;
    ylabel("Desplazamiento [in]");
    xlabel("Tiempo [s]");
    grid on;

    yyaxis right;
    hold on;
    plot(t,Fk,"LineWidth",1.5)
    % scatter(t,P,"filled");
    hold off;
    ylabel("Fuerza Elastica [kips]");
    xlabel("Tiempo [s]");
    grid on;
    legend("X(t)","Fk(t)");

    % mostar una tabla con los primeros 10 n valores de P, x ,dx y Fk
    n = 0:1:10;
    T = table(n(1:10)',P(1:10)',xt(1:10)',dx(1:10)',Fk(1:10)','VariableNames',{'n','Pn','x(t)','v(t)','Fk(t)'});
    disp(T);

end

function P=carga(t) % Carga externa en el tiempo
    t_app = 0:0.12:0.72;
    P_app = [0,1,4,9,9,6,0];
    P = interp1(t_app,P_app,t,"linear",0);
end

function [x,dx] = PaP_dif_central(T,t,dt,P,m,c,k)
    if dt/T<=1/pi
        %Inicializo vectores de desplazamiento, velocidad
        x=zeros(1,length(t));
        dx=zeros(1,length(t));

        %Establezco condiciones iniciales
        x(1)=0;
        dx(1)=0;

        %Aplico ecuaciones de Newmark Beta
        for i=2:length(t) %ojo, calculo del instante 2 en adelante
            x(i)=x(i-1)+dt*dx(i-1)+(dt^2)*(P(i-1)-c*dx(i-1)-k*x(i-1))/(2*m);
            dx(i)=2*(x(i)-x(i-1))/dt-dx(i-1);
        end
    else
        disp('Error: dt/T debe ser menor o igual a 1/pi');
        disp(['dt/T: ',num2str(dt/T)]);
        disp(['1/pi: ',num2str(1/pi)]);
    end
    
end