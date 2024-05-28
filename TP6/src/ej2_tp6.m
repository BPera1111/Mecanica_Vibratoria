function ej2_tp6;clc;close all;
    % Parametros del sistema con un grado de libertad
    m = 250; % Masa [kg]
    k = 4000e4; % Rigidez [N/m]
    zitta = 0.05; % Coeficiente de amortiguamiento   
    wn = sqrt(k/m); % Frecuencia natural [rad/s]

    wd = wn*sqrt(1-zitta^2); % Frecuencia de amortiguamiento [rad/s]

    dt = 0.001; % Paso de tiempo [s]
    % tf = 5; % Tiempo final [s]
    t = 0:dt:0.3; % Vector de tiempo [s]

    % Parametros de la carga externa
    Po = 430e3; % Amplitud de carga externa [kN/m]
    Tp = 0.05; % Periodo de carga externa [s]
    P = carga(t,Tp,Po); % Carga externa en el tiempo

    xt = Duhamel(t,dt,P,wd,zitta,m,k); % Calculo de la integral de Duhamel

    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo');
    yyaxis left;
    hold on;
    plot(t,xt,"LineWidth",1.5)
    % scatter(t,xt,"filled");
    hold off;
    ylabel("x(t) [m]");

    yyaxis right;
    hold on;
    plot(t,P,"LineWidth",1.5)
    % scatter(t,P,"filled");
    hold off;
    ylabel("P(t) [N]");
    xlabel("Tiempo [s]");
    grid on;
    legend("Desplazamiento","Carga");
end


function P = carga(t,tp,Po) % Carga externa en el tiempo
    P = zeros(1,length(t));
    for i = 1:length(t)
        % t_mod = mod(t(i),tp);
        P(i) = (t(i)<=tp/2).*(t(i)*Po/(tp/2)) + (t(i)>tp/2 && t(i)<=tp).*(t(i)*(-Po)/(tp/2)+2*Po);
    end
end

function xt=Duhamel(t,dt,P,wd,zitta,m,k) % Integral de Duhamel para sistemas subamortiguados
    e=exp(1)^(-zitta*wd*dt);
    A = zeros(1,length(t));
    B = zeros(1,length(t));
    yc = P.*cos(wd*t);
    ys = P.*sin(wd*t);
    for i=2:length(t)
        for j=1:i
            A(i) = A(i-1)*e+ dt/(2*m*wd) * (yc(i-1)*e + yc(i));
            B(i) = B(i-1)*e+ dt/(2*m*wd) * (ys(i-1)*e + ys(i));
        end
    end
    xt = A.*sin(wd*t) - B.*cos(wd*t);

    %fuerza del resorte
    Fk = (-k)*xt;

    % tabla con las primeras 10 filas de P, yc, ys, A , B , xt y Fk
    n = 0:1:10;
    T = table(n(1:10)',P(1:10)',yc(1:10)',ys(1:10)',A(1:10)',B(1:10)',xt(1:10)',Fk(1:10)','VariableNames',{'n','Pn','yc(t)','ys(t)','A(t)','B(t)','x(t)','Fk(t)'});
    disp(T); 
end