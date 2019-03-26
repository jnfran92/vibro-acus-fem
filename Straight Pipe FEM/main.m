
% Tarea #3
% ACUS 360
% Modelamiento de los modos acusticos de un tubo por elementos finitos
% Autor: Juan Chango

clc;
tic                                 % inicio de medicion de tiempo


% Variables de Entrada ------------------------------------------------
% Solo se requiere modificar las siguientes variables

n_elements = 64;                   % numero de elementos equidistantes en la viga
n_dofs = 1;                         % grados de libertad
n_nodes = n_elements + 1;           % numero de nodos equidistantes en la viga

cons = [];                          % grados de libertad iguales a cero
                                    % ingresar como vector

rho = 1.21;                         % Densidad
c = 343;                            % Velocidad del Somido

L = 0.5;                            % Largo de la viga
A = 0.05;                           % Seccion de la viga

%-----------------------------------------------------------------------


% Calculo de las matrices de rgidez e inercia
% del elemento individual:

a1 =  L/(n_nodes-1)/2;
n_nodes_element = 2;

Q = zeros(n_nodes_element * n_dofs); 
H = zeros(n_nodes_element * n_dofs);


q_cons = (A * a1) / c^2;
Q = q_cons * [  2/3       1/3; ...
                1/3       2/3 ] ;



h_cons = A / a1;
H = h_cons * [  1/2      -1/2  ; ...
               -1/2       1/2    ] ;


            
% Creacion de Matrices globales de inercia y rigidez  

Qg = zeros(n_nodes * n_dofs);       % Matriz global de inercia
Hg = zeros(n_nodes * n_dofs);       % Matriz global de rigidez
            
p = n_nodes_element * n_dofs;

% Ensamblaje de la matriz global usango matriz auxiliar Ae
for e=1:n_nodes - 1
    Ae =  [ zeros(p,1*(e-1)),...    
            eye(p),...
            zeros(p, n_nodes * n_dofs - p - 1*(e-1))];

    Qg = Qg + Ae' * Q * Ae;
    Hg = Hg + Ae' * H * Ae;
end


% Eliminacion de grados de libertad iguales a zero (vector 'cons')
Qg(:,cons) = [];
Qg(cons,:) = [];

Hg(:,cons) = [];
Hg(cons,:) = [];


% Calculo de los vectores y valores propios
[V,LAMBDAI] = eig(Hg,Qg);
LAMBDA = LAMBDAI*ones(length(LAMBDAI),1);


% Resultados
disp('*****************************************');
disp(' ');
disp('Modelamiento de los Modos Acusticos')
disp('de un Tubo por Elementos Finitos')
disp(' ');
disp('*****************************************');
disp('Caracteristicas del tubo');
disp(strcat('Largo(m): ',num2str(L) ))
disp(strcat('Seccion transversal(m^2): ',num2str(A) ))

disp('*****************************************');
disp('Detalles del modelamiento');
disp(strcat('Numero total de nodos: ',num2str(n_nodes) ))
disp(strcat('Numero total de elementos: ',num2str(n_elements) ))

disp('*****************************************');
disp('Resultados');
% Frecuencias naturales del sistema
Freqs = round(( LAMBDA.^(0.5) )./(2*pi));
% Se imprimen en consola las frecuencias
disp('5 Primeras Frecuencias naturales')
disp(Freqs(1:5))
disp('Tiempo de modelamiento')
toc


% Grafico de las 4 frecuencias mas bajas

ButtonName = questdlg('Desea graficar los resultados?');

if strcmp(ButtonName, 'Yes')
    [X,Y] = meshgrid(0:n_elements,0:1);
    
    for i=1:5
        subplot(5,1,i);
        hold on
        
        index = i +1 ;
        n_freq = i +1;
        
        
        C =  [V(:,n_freq)';V(:,n_freq)'];
        surf(X,Y,X*0,C)
        view(2)
        axis([0 n_elements -1 1])
        colormap('jet')
        
        
        if i == 1
            title_str = { strcat('SIMULACION PARA: ',...
                num2str(n_elements),...
                ' ELEMENTOS');...
                strcat('Modo #',num2str(index-1),...
                '-Frecuencia:',num2str(Freqs(index)), ' Hz')};
        else
            title_str = strcat('Modo #',num2str(index-1),...
                '-Frecuencia:',num2str(Freqs(index)), ' Hz');
        end
        
        
        title(title_str)
        axis off
    end
    
end



ButtonName2 = questdlg('Desea generar una animacion a partir de los resultados(Atencion: El proceso puede tomar mucho tiempo)?');

% Crear Animacion 
n_frames = 80;
n_loop = 2;
counter = 1;


if strcmp(ButtonName2, 'Yes')
    
    figure
    for i=1:4
        n_freq = i +1;
        index = i +1 ;
        
        for j= 1:n_loop
            for k = 0:n_frames
                
                f = cos((k/n_frames)*(2)*pi);
                C =  [(f)*V(:,n_freq)';(f)*V(:,n_freq)'];
                surf(X,Y,X*0,C)
                caxis([-700, 700]);
                axis off;
                cb = colorbar('south');
                title(cb,'Presion Sonora')
                
                title_str = { strcat('Simulacion con: ',...
                    num2str(n_elements),...
                    ' Elementos');...
                    strcat('Modo #',num2str(index-1),...
                    '-Frecuencia:',num2str(Freqs(index)), ' Hz')
                    };
                
                dim = [.35 .6 .3 .3];
                an = annotation('textbox',dim,'String',title_str,'FitBoxToText','on',...
                    'BackgroundColor', 'white');
                an.FontSize = 14;
                
                
                view(2)
                axis([0 n_elements -1 2])
                colormap('jet')
                
                
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
end




