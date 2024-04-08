function ej4; clc; close all; clear all;

    m = 1361;
    k = 2.688E5;
    c = 3.81E3;

    w_n=sqrt(k/m);
    zeta=c/(2*m*w_n);
    w_d=w_n*sqrt(1-zeta^2);

    t = 0:0.01:4;
    x_0 = 0;
    v_0 = 0.01;

    A = x_0;
    B = (v_0 + zeta*w_n*x_0)/w_d;

    % solución analítica

    x = A*exp(-zeta*w_n*t).*cos(w_d*t) + B*exp(-zeta*w_n*t).*sin(w_d*t);

    % solución numérica

    function f = sys(t,x)
    %sys ecuación diferencial de un sistema de segundo orden
    
        m = 1361;
        k = 2.688E5;
        c = 3.81E3;

        f = [x(2); -c/m*x(2) - k/m*x(1)];

    end

    [t,y] = ode45(@sys,t,[x_0;v_0]);

    figure;
    hold on;
    plot(t,y(:,1),'b');
    plot(t,x,'ro');
    legend('Solución numérica','Solución analítica');
    xlabel('t');
    ylabel('x(t)');
    title('Respuesta al impulso de un sistema de segundo orden');
    grid on;



end