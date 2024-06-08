function ej4_tp5 ;clc;close all;
    
    % Parametros de la carga externa

    Tp = 2*pi; % Periodo de carga externa [s]
    wp = 2*pi/Tp; % Frecuencia de carga externa [rad/s]
    Po = 1; % Amplitud de carga externa [N]
    dt = 0.0001; % Paso de tiempo [s]

    t = 0:dt:Tp; % Vector de tiempo [s]
    tg = 0:dt:5*Tp; % Vector de tiempo extendido [s]

    Pt = carga(t,Tp,Po); % Carga externa [N]
    Ptg = carga(tg,Tp,Po); % Carga externa extendida [N]

    % Mostramos la carga externa orginal en el tiempo
    figure;
    subplot(2,2,1);
    plot(tg,Ptg);
    title("Carga externa en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("Carga [N]");
    grid on;

    % Carga externa con Fourier

    n = 10; % Cantidad de armonicos

    % Calculo de los coeficientes de Fourier
    a0 = 2/Tp*trapz(t,Pt); % Coeficiente de Fourier a0

    A = zeros(length(t),n); % Coeficientes de Fourier A
    B = zeros(length(t),n); % Coeficientes de Fourier B

    an = zeros(n,1); % Coeficientes de Fourier an
    bn = zeros(n,1); % Coeficientes de Fourier bn

    % Calculo de los coeficientes de Fourier
    for j = 1:n
        A(:,j) = Pt.*cos(j*wp*t); % Coeficientes de Fourier A
        B(:,j) = Pt.*sin(j*wp*t); % Coeficientes de Fourier B
        an(j) = 2/Tp*trapz(t,A(:,j)); % Coeficientes de Fourier an
        bn(j) = 2/Tp*trapz(t,B(:,j)); % Coeficientes de Fourier bn
    end
    a0
    an
    bn

    % Calculo de la carga externa con Fourier
    % Armonicos
    armonicosg = armo(tg,wp,a0,an,bn,n);
    % armonicos = armo(t,wp,a0,an,bn,n);

    Pt_fourierg2 = carga_f2(n,armonicosg);
    % Pt_fourier2 = carga_f2(n,armonicos);

    % Mostramos los armonicos en el tiempo
    subplot(2,2,2);
    plot(tg,armonicosg);
    title("Armonicos en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("Carga [N]");
    grid on;

    % Mostramos la carga externa con Fourier en el tiempo
    subplot(2,2,3);
    plot(tg,Pt_fourierg2);
    title("Carga externa con Fourier en el tiempo");
    xlabel("Tiempo [s]");
    ylabel("Carga [N]");
    grid on;

    % Respuesta en el tiempo

    % Parametros del sistema
    m = 100; % Masa [kg]
    k = 1000; % Rigidez [N/m]
    w = sqrt(k/m); % Frecuencia natural [rad/s]
    zitta = 0.1; % Coeficiente de amortiguamiento
    betta = wp/w; % Relacion de frecuencias
    


    % Calculo de la respuesta en el tiempo
    xt = respuesta(tg,a0,an,bn,n,wp,k,zitta,betta);
    
    

    % Mostramos la respuesta en el tiempo
    subplot(2,2,4);
    [hAx,hLine1,hLine2] = plotyy(tg,Pt_fourierg2,tg,xt);
    set(hLine1,'LineStyle','-','LineWidth',1.5,'Color','r');
    set(hLine2,'LineStyle','-','LineWidth',1.5,'Color','b');
    title("Respuesta en el tiempo");
    xlabel("Tiempo [s]");
    ylabel(hAx(1),"Carga [N]");
    ylabel(hAx(2),"Desplazamiento [m]");
    grid on;



end

function Pt =carga(t,tp,Po) % Calculo de la carga externa en el tiempo
    Pt = zeros(1,length(t));
    for i = 1:length(t)
        t_mod = mod(t(i),tp);
        Pt(i) = (t_mod>=0 && t_mod<tp/2).*(4/tp*t_mod-Po) + (t_mod>=tp/2 && t_mod<tp).*(Po-4/tp*(t_mod-tp/2));
    end

end


function Pt_fourier = carga_f2(n,armonicos) % Calculo de la carga externa con Fourier
    Pt_fourier = zeros(1,length(armonicos)); % Definimos el tamaño de la carga externa con Fourier [N] en funcion de los armonicos
    for j = 1:n
        Pt_fourier = Pt_fourier + armonicos(j,:); % Realizamos la suma de los armonicos
    end
end

function arm = armo(t,wp,a0,an,bn,n) % Calculo de los armonicos
    arm = zeros(n,length(t)); % Definimos el tamaño de los armonicos en funcion de la cantidad armonicos y el tiempo
    for i = 1:length(t)
        for j = 1:n
            arm(j,i) = a0/2 + an(j)*cos(j*wp*t(i)) + bn(j)*sin(j*wp*t(i)); % Se calcula el valor de cada armonico en cada instante de tiempo
        end
    end
end

function xt = respuesta(tg,a0,an,bn,n,wp,k,zitta,betta) % Calculo de la respuesta en el tiempo
    phi = zeros(n,1); % Fase de cada armonico
    H = zeros(n,1); % Funcion de transferencia de cada armonico

    for j = 1:n
        phi(j) = atan((2*zitta*betta*j)/(1-betta^2*j^2)); % Fase de cada armonico
        H(j) = 1/(k*sqrt((1-j^2*betta^2)^2+(2*zitta*j*betta)^2)); % Funcion de transferencia de cada armonico
    end

    xt = zeros(1,length(tg)); % vector que almacena la respuesta en el tiempo 

    for i = 1:length(tg)
        for j = 1:n
            xt(i) = xt(i) + a0/(2*k)+ an(j)*H(j)*cos(j*wp*tg(i)-phi(j)) + bn(j)*H(j)*sin(j*wp*tg(i)-phi(j)); % Se calcula el valor de cada armonico en cada instante de tiempo
        end
    end
end