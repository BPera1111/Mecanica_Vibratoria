function MultiplesGradosGeneral
 clc,clear
 carga = 'L';
 dt = 0.01;
 t = 0:dt:2;
 %Amortiguamiento
 c1 = 0;
 c2 = 0;
 c3 = 0;
 %Masas
 m1 = 2;
 m2 = 2;
 m3 = 2;
 %Rigidez
 k1 = 600;
 k2 = 1200;
 k3 = 2400;

 %Condiciones iniciales
x0 = [0.3; -0.8 ; 0.3]; %Desplazamientos iniciales
xd0 = [0; 0; 0] ; %Velocidades iniciales

 m= diag([m1 m2 m3]); %Matriz de masas
 k = [k1, -k1, 0;
        -k1, k1+k2, -k2;
          0, -k2, k2+k3]; %Matriz de rigidez
c = zeros(size(m)); %Constantes de amort

I = eye (size(m));
[V, lam] = eig(k,m); %Calculo de autovalores y autovectores
ome = diag(lam).^0.5; % Los omega o fecuencias naturales de cada modo
ttiemp = 2*pi/ome(1)
idx = find (t == ttiemp)

M  = V' * m * V;  %Matriz de masas modales
V
M
I
V'
m
x0
Y0  = M\I*V'*m*x0  %Desplazamientos en coordenadas modales
Yd0  = M\I*V'*m*xd0; %Velocidades en coordenadas modales
M\I*V'*m
%Acá lo que aplico es la ecuacion general de la respuesta de un sistema libre. x(t) = a * cos... + b * sin ...
%Donde a = x0 y b = xd0/wn (sin amortiguamiento)
%Primer término de la solucion en coordenadas modales
 zitta1 = c1/(2*ome(1)*m1) ;
 zitta2 = c2/(2*ome(2)*m2) ;
 zitta3 = c3/(2*ome(3)*m3) ;

 if (carga == 'L')

   Y1  = [Y0(1)*cos(sqrt(1-zitta1^2)*ome(1)*t);
            Y0(2)*cos(sqrt(1-zitta2^2)*ome(2)*t);
            Y0(3)*cos(sqrt(1-zitta3^2)*ome(3)*t)];
    %Segundo termino de la solucion en coordenadas modales
    C21 = (Yd0(1)+(zitta1*ome(1)*Y0(1)))/(sqrt(1-zitta1^2)*ome(1));
    C22 = (Yd0(2)+(zitta2*ome(2)*Y0(2)))/(sqrt(1-zitta2^2)*ome(2));
    C23 = (Yd0(3)+(zitta3*ome(3)*Y0(3)))/(sqrt(1-zitta3^2)*ome(3));

    Y2  = [C21*sin(sqrt(1-zitta1^2)*ome(1)*t);
                C22*sin(sqrt(1-zitta2^2)*ome(2)*t);
                C23*sin(sqrt(1-zitta3^2)*ome(3)*t)];
    e = exp(1)
    ee = [e.^(-zitta1*ome(1)*t);
              e.^(-zitta2*ome(2)*t);
              e.^(-zitta3*ome(3)*t)];
    Y = (Y1 + Y2).*ee; %Desplazamientos en coords modales
    X = V*Y; %Lo pasamos a coordenadas geometricas generalizadas

    subplot(3,1,1)
    plot(t, X(1,:), 'b');

    subplot(3,1,2)
    plot(t, X(2,:), 'r');

    subplot(3,1,3)
    plot(t, X(3,:), 'g');
 end

end
