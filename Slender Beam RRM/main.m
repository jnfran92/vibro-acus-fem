
clc;
% Variables de Entrada del Sistema
n_deg_free = 8; % grados de libertad

rho = 2700; % densidad
L = 0.5; % Largo de la viga
E =  7.1E10; % Modulo de Young

b = 0.02; % ancho
h = 0.005; % grosor

A = b * h; % Secci?n de la viga
Iz = (b*h^3)/12; % Segundo Momento polar Rectangulo con eje pasando por el centro de gravedad

% Inicilizando matrices de masa e inercia
M = zeros(n_deg_free); 
K = zeros(n_deg_free);

% C?lculo de las matrices de M y K
for k=1:n_deg_free
    for s=1:n_deg_free
        
        m1 = (L^(k+s+5))/(k+s+3);
        m2 = 2*((L^(k+s+5))/(k+s+4));
        m3 = ((L^(k+s+5))/(k+s+5));
        
        M(k,s) = rho*A*( m1-m2+m3 );
        
        
        k1 = ( ( (k+1)*k*(s+1)*s*L^(k+s+1) ) / (k+s-1));
        k2 = ( ( (k+2)*(k+1)*(s+1)*s*L^(k+s+1) ) / (k+s));
        k3 = ( ( (k+1)*k*(s+2)*(s+1)*L^(k+s+1) ) / (k+s));
        k4 = ( ( (k+2)*(k+1)*(s+2)*(s+1)*L^(k+s+1) ) / (k+s+1));
        
        K(k,s) = E*Iz*( k1-k2-k3+k4 );
    end
    
end

% C?lculos de los vectores y valores propios
[A,LAMBDAI] = eig(K,M);
LAMBDA = LAMBDAI*ones(n_deg_free,1);


% Frecuencias naturales del sistema
Freqs = round(( LAMBDA.^(0.5) )./(2*pi));
% Se imprimen en consola las frecuencias
disp('Frecuencias naturales')
disp(Freqs)

% Caculo de los resultados
n_points = 1000;
x = 0:L/(n_points-1):L;
d = zeros(n_deg_free,n_points);

% Calculo de las funciones espaciales de la ecuacion de Rayleigh-Ritz
for i=1:n_deg_free 
    d(i,:) = (x.^(i+1)).*(L-x);
end

% Matriz de resultados
Y = (d'*A)';


% Graficos de los resultados
legs = {};
for i=1:4
   plot(x, Y(i,:),'LineWidth', 1.6);
   legs{i} =  ['Modo #:', num2str(i),' Frecuencia:', num2str(Freqs(i)),'Hz'];
   hold on
end
grid on
title(['Modos Naturales de Vibracion: ',num2str(n_deg_free), ' grados de libertad']);
lgd = legend(legs);
lgd.FontSize = 14;


% Animacion 
vidObj = VideoWriter('animacion.mp4','MPEG-4');
open(vidObj);
zs = zeros(1,n_points);
figure 
loops = 10;
for i=1:n_deg_free
    F(loops) = struct('cdata',[],'colormap',[]);
    for j = 1:loops
        plot(x, Y(i,:),'LineWidth', 1.6, 'Color', 'r')
        hold on
        plot(x,zs,'LineWidth', 1.9, 'Color', 'b')
        title(['Modo Natural de Vibracion #',num2str(i)]);
        lgd = legend(['Modo #',num2str(i),'  Frecuencia:', num2str(Freqs(i)),'Hz']);
        lgd.FontSize = 14;
        axis off
        ylim([min(min(Y)) max(max(Y))])
        drawnow
        F(j) = getframe;
        writeVideo(vidObj,F(j));
    end
    hold off
end

% Guardar animacion
close(vidObj);


