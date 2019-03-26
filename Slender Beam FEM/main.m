clear all
clc;
% Elemento finito viga delgada
% Autor: Juan Chango
tic


n_elements = 4;                     % numero de elementos equidistantes en la viga

% Variables de Entrada
n_dofs = 2;                         % grados de libertad
n_nodes = n_elements + 1;           % numero de nodos equidistantes en la viga

cons = [1,2,n_nodes * n_dofs-1];    % grados de libertad iguales a cero
 

rho = 2700;                         % Densidad
E =  7.1E10;                        % Modulo de Young

L = 0.5;                            % Largo de la viga
b = 0.02;                           % Ancho
h = 0.005;                          % Grosor


% Calculo de propiedads geom?tricas de la viga

A = b * h;                          % Seccion de la viga
Iz = (b*h^3)/12;                    % Segundo Momento polar Rectangulo con 
                                    %eje pasando por el centro de gravedad

% Calculo de las matrices de masa e inercia 
% del elemento individual:

a =  L/(n_nodes-1) / 2;
n_nodes_element = 2;

M = zeros(n_nodes_element * n_dofs); 
K = zeros(n_nodes_element * n_dofs);


m_cons = (rho * A * a) / (105);
M = m_cons * [  78      22*a    27      -13*a; ...
                22*a    8*a^2   13*a    -6*a^2;...
                27      13*a    78      -22*a;...
                -13*a   -6*a^2 -22*a     8*a^2 ] ;



k_cons = (E * Iz) / (2 * a^3);
K = k_cons * [  3       3*a     -3      3*a; ...
                3*a     4*a^2   -3*a    2*a^2;...
                -3      -3*a     3      -3*a;...
                3*a     2*a^2   -3*a    4*a^2 ] ;


            
% Creaci?n de Matrices globales de masa y rigidez  

Mg = zeros(n_nodes * n_dofs);       % Matriz global de inercia
Kg = zeros(n_nodes * n_dofs);       % Matriz global de rigidez
            
p = n_nodes_element * n_dofs;

% Ensamblaje de la matriz global usango matriz auxiliar Ae
for e=1:n_nodes - 1
    Ae =  [ zeros(p,2*(e-1)),...    
            eye(p),...
            zeros(p, n_nodes * n_dofs - p - 2*(e-1))];

    Mg = Mg + Ae' * M * Ae;
    Kg = Kg + Ae' * K * Ae;
end


% Eliminacion de grados de libertad iguales a zero (vector 'cons')
Mg(:,cons) = [];
Mg(cons,:) = [];

Kg(:,cons) = [];
Kg(cons,:) = [];


% Calculo de los vectores y valores propios
[A,LAMBDAI] = eig(Kg,Mg);
LAMBDA = LAMBDAI*ones(length(LAMBDAI),1);


% Resultados

% Frecuencias naturales del sistema
Freqs = round(( LAMBDA.^(0.5) )./(2*pi));
% Se imprimen en consola las frecuencias
disp('Frecuencias naturales')
disp(Freqs)

toc


%Granficando los modos y grabando animaciones
idx = 1:n_nodes * n_dofs;
idx(cons) = [];
R = zeros(n_nodes * n_dofs,length(idx)); 
R(idx,:) = A;

% Gr?fico de las 4 frecuencias m?s basjas 
for i=1:4
subplot(4,1,i);
hold on
index = i;
aux = R(1:2:end,index);
plot(aux)
plot(aux,'*')
plot(zeros(1,length(aux)))
if i == 1
    title_str = { strcat('Simulacion para: ',...
        num2str(n_elements),...
        ' elementos');...
        strcat(num2str(Freqs(index)), ' Hz')};
else
    title_str = strcat(num2str(Freqs(index)), ' Hz');
end
title(title_str)
axis off
end



%Animaci?n------------------------------------------
n_frames = 20;
n_loop = 2;
counter = 1;

for i=1:4
    index = i
    for j= 1:n_loop
        for k = -n_frames:n_frames
            
            aux = (k/n_frames)*R(1:2:end,index);
            plot(aux)
            hold on
            plot(aux,'*')
            plot(zeros(1,length(aux)))
            title_str = { strcat('Simulacion para: ',...
                num2str(n_elements),...
                ' elementos  ',...
                   num2str(Freqs(index)), ' Hz')};
            legend(title_str)
            axis([0 n_elements+2 -6 6])
            T(counter) = getframe;
            counter = counter +1;
            
            hold off
        end
        
        for k = -n_frames:n_frames
            
            aux = (-k/n_frames)*R(1:2:end,index);
            plot(aux)
            hold on
            plot(aux,'*')
            plot(zeros(1,length(aux)))
            
             title_str = { strcat('Simulacion para: ',...
                num2str(n_elements),...
                ' elementos ',...
                   num2str(Freqs(index)), ' Hz')};
            legend(title_str)
            axis([0 n_elements+2 -6 6])
            T(counter) = getframe;
            counter = counter +1;
            
            hold off
        end
        
    end
end

vidObj = VideoWriter(strcat('anim_',...
    num2str(n_elements),'elementos.mp4'),'MPEG-4');
open(vidObj);
writeVideo(vidObj,T);
close(vidObj);



