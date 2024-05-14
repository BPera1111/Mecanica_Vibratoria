function ej1_tp6 ;clc;close all;

    % Parametros del sistema con un grado de libertad
    m = 0.2; % Masa [kips.sec^2/in]
    k = 8; % Rigidez [kips/in]
    c = 0.4; % Amortiguamiento [kips.sec/in]
    wn = sqrt(k/m); % Frecuencia natural [rad/s]
    zitta = c/(2*sqrt(m*k)); % Coeficiente de amortiguamiento
    wd = wn*sqrt(1-zitta^2); % Frecuencia de amortiguamiento [rad/s]

    dt = 0.12; % Paso de tiempo [s]
    tf = 5; % Tiempo final [s]
    t = 0:dt:tf; % Vector de tiempo [s]

    % Parametros de la carga externa
    P = carga(t); % Carga externa en el tiempo

    % Calculo de la integral de Duhamel
    xt = Duhamel(t,dt,P,wd,zitta,m,k);

    % Calculo de la respuesta en el tiempo
    % xt = A.*sin(wd*t) - B.*cos(wd*t);

    figure('Name', 'Respuesta en el tiempo superpuesta');
    plot(t,P,'r','LineWidth',1.5);
    hold on;
    plot(t,xt,'b','LineWidth',1.5);
    scatter(t,P,'filled','r');
    scatter(t,xt,'filled','b');
    hold off;
    title("Respuesta en el tiempo superpuesta");
    xlabel("Tiempo [s]");
    ylabel("Carga [kips] / Desplazamiento [in]");
    grid on;

    


    % Mostramos la respuesta en el tiempo
    figure('Name', 'Respuesta en el tiempo');
    [hAx,hLine1,hLine2] = plotyy(t,P, t, xt);
    set(hLine1,'LineStyle','-','LineWidth',1.5,'Color','r');
    set(hLine2,'LineStyle','-','LineWidth',1.5,'Color','b');
    hold(hAx(1),'on');
    scatter(hAx(1),t,P,'filled','r');
    hold(hAx(1),'off');
    hold(hAx(2),'on');
    scatter(hAx(2),t,xt,'filled','b');
    hold(hAx(2),'off');
    title("Respuesta en el tiempo");
    xlabel("Tiempo [s]");
    ylabel(hAx(1),"Carga [kips]");
    ylabel(hAx(2),"Desplazamiento [in]");
    grid on;

end

function P=carga(t) % Carga externa en el tiempo
    t_app = 0:0.12:0.72;
    P_app = [0,1,4,9,9,6,0];
    P = interp1(t_app,P_app,t,"linear",0);
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
