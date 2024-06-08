function ej3_tp5; clc;close all;
    % Parametros del sistema con un grado de libertad
    m = 1; % Masa [kg]
    k = 1; % Rigidez [N/m]
    wn = sqrt(k/m); % Frecuencia natural [rad/s]
    T = 2*pi/wn; % Periodo [s]
    

    % Parametros de la carga externa
    Po = 1; % Amplitud de carga externa [N]
    Tp = T*4/3; % Periodo de carga externa [s]
    wp = 2*pi/Tp; % Frecuencia de carga externa [rad/s]
    dt = 30*pi/(wp*180); % Paso de tiempo [s]
    tg = 0:dt:Tp*2; % Vector de tiempo [s]

    zitta = 0; % Coeficiente de amortiguamiento
    betta = wp/wn; % Relacion de frecuencias

    n = 4; % Cantidad de armonicos

    %calculo de los coeficientes de Fourier

    bn = zeros(n,1); % Coeficientes de Fourier bn
    for j = 1:n
        bn(j) = -(2*Po)/(pi*j)*(-1)^j; % Coeficiente de Fourier bn
    end

    % Pt = carga(tg,n,Po,wp); % Carga externa en el tiempo
    armonicos = armo(tg,wp,bn,n); % Armonicos
    Pt_fourier = carga_f2(n,armonicos); % Carga externa con Fourier
    Pt = carga(tg,Tp,Po); % Carga externa en el tiempo
    xt = respuesta(tg,bn,n,wp,k,zitta,betta); % Respuesta en el tiempo

    Pt_fourier(1:13)
    % Mostramos la carga externa en el tiempo
    figure;
    subplot(2,2,1);
    plot(tg,Pt);
    hold on
    scatter(tg,Pt,'filled');
    hold off
    title("Carga externa en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("P(t) [N]");
    grid on;

    % Mostramos los armonicos en el tiempo
    subplot(2,2,2);
    plot(tg,armonicos);
    title("Armonicos en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("Carga [N]");
    grid on;

    % Mostramos la carga externa con Fourier en el tiempo
    subplot(2,2,3);
    plot(tg,Pt_fourier);
    hold on
    scatter(tg,Pt_fourier,'filled');
    hold off
    title("Carga externa con Fourier en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("Carga [N]");
    grid on;

    % Mostramos la respuesta en el tiempo
    subplot(2,2,4);
    plot(tg,xt);
    hold on
    scatter(tg,xt,'filled');
    hold off
    title("Respuesta en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("x(t) [m]");
    grid on;

end

function Pt =carga(t,tp,Po) % Calculo de la carga externa en el tiempo
    Pt = zeros(1,length(t));
    for i = 1:length(t)
        t_mod = mod(t(i),tp); % mod() devuelve el resto de la division de t(i) por tp
        Pt(i) = (t_mod>=0 && t_mod<tp/2).*(2/tp*t_mod) + (t_mod>=tp/2 && t_mod<tp).*(2/tp*t_mod-2*Po);
    end

end


function Pt_fourier = carga_f2(n,armonicos) % Calculo de la carga externa con Fourier
    Pt_fourier = zeros(1,length(armonicos)); % Definimos el tamaño de la carga externa con Fourier [N] en funcion de los armonicos
    for j = 1:n
        Pt_fourier = Pt_fourier + armonicos(j,:); % Realizamos la suma de los armonicos
    end
end

function arm = armo(t,wp,bn,n) % Calculo de los armonicos
    arm = zeros(n,length(t)); % Definimos el tamaño de los armonicos en funcion de la cantidad armonicos y el tiempo
    for i = 1:length(t)
        for j = 1:n
            arm(j,i) =  bn(j)*sin(j*wp*t(i)); % Se calcula el valor de cada armonico en cada instante de tiempo
        end
    end
end

function xt = respuesta(tg,bn,n,wp,k,zitta,betta) % Calculo de la respuesta en el tiempo
    phi = zeros(n,1); % Fase de cada armonico
    H = zeros(n,1); % Funcion de transferencia de cada armonico

    for j = 1:n
        phi(j) = atan((2*zitta*betta*j)/(1-betta^2*j^2)); % Fase de cada armonico
        H(j) = 1/(k*sqrt((1-j^2*betta^2)^2+(2*zitta*j*betta)^2)); % Funcion de transferencia de cada armonico
    end

    xt = zeros(1,length(tg)); % vector que almacena la respuesta en el tiempo 

    for i = 1:length(tg)
        for j = 1:n
            xt(i) = xt(i) + bn(j)*H(j)*sin(j*wp*tg(i)-phi(j)); % Se calcula el valor de cada armonico en cada instante de tiempo
        end
    end
end