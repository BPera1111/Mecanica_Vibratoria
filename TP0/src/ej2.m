function ej2

t=[0:2*pi/50:2*pi];
x=2*cos(t);
y=5*sin(t);

figure;

subplot(2,2,1);
title('Gráficos de x(t) y y(t)');
xlabel('t'); ylabel('x(t), y(t)');
hold on;
plot(t,x);
plot(t,y);
hold off;

subplot(2,2,2);
title('Gráfico de x(t) vs y(t)');
xlabel('x(t)'); ylabel('y(t)');
plot(x,y);

subplot(2,2,3);
title('Gráfico de x(t) vs y(t) con ejes iguales');
xlabel('x(t)'); ylabel('y(t)');
hold on;
axis equal; %para que los ejes tengan la misma escala
plot(x,y)
hold off;

subplot(2,2,4);
paso = input('Ingrese el paso de t para la animación ej(2*pi/30): ');
t2=[0:paso:2*pi];
title('Animación');
axis([-5 5 -5 5])
hold on; grid on;
xlabel('X(t)'); ylabel('Y(t)');
for i=1:length(t2)-1
    plot(x(i:i+1),y(i:i+1),'LineWidth',2);
    pause(0.1);
end
hold off;

end
    



