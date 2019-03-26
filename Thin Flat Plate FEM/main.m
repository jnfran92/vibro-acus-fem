
% Tarea #4
% ACUS 360
% Modelamiento de los Modos de Flexion de una Placa Fina
% Autor: Juan Chango

clc
clear all                       % Borrar todas las variables

% Variables de Entrada ------------------------------------------------
% Solo se requiere modificar las siguientes variables

nd = 33;                        % Grilla de tama?o nd x nd, 
                                % nd:numero de divisiones
 
Lx = 0.414;                     % Largo
Ly = 0.314;                     % Ancho
h = 0.001;                      % Espesor

densidad_rho = 2700;            % Densidad del material Kg / m^3
modulo_young = 7.1e10;          % Modulo de Young
coef_poisson =  0.33;           % Coeficiente de Poisson

%-----------------------------------------------------------------------

%   1. Generaci?n de ecuaciones con variables simb?licas para
%   encontrar las matrices de inercia y rigidez, las ecuaciones de 
%   referencia son de: F. J. Fahy and P. Gardonio, Sound and structural 
%   vibration: radiation, transmission and response. Elsevier, 2007.

syms E n                        % Variables: xi y eta del polinomio
syms a b                        % a y b largo y ancho del elemento
syms Em v                       % Modulo de Young y Coef. Poisson

% Polinomio vector y transpuesto ref. Eq. (8.45a)
p_t = [
    1;...
    E;...
    n;...
    E^2;...
    E*n;...
    n^2;...
    E^3;...
    (E^2)*n;...
    E*(n^2);...
    n^3;...
    (E^3)*n;...
    E*(n^3)
    ];

p = [
    1,...
    E,...
    n,...
    E^2,...
    E*n,...
    n^2,...
    E^3,...
    (E^2)*n,...
    E*(n^2),...
    n^3,...
    (E^3)*n,...
    E*(n^3)
    ];

% Multiplicacion de [p]^T*[p] e Integraci?n para obtener matriz de masa
% ref. Eq. (8.51a)
pp_m = p_t*p;
pp_m_int = int(int(pp_m,E,[-1 1]),n,[-1 1]);


% Generaci?n de la matriz de rigidez
% Em_p: E' y Em: Modulo de Young ref. Eq. (8.41b)
Em_p = Em/(1 - v^2);
G = Em/(2* (1 + v));

D = [ 
    Em_p     Em_p*v    0   ;...
    Em_p*v    Em_p     0   ;...
    0         0        G
    ];

% Matriz de polinomios(segunda derivada) ref. Eq.(8.52)
P_t = [
    (1/a^2) * diff(diff(p_t,E),E),...
    (1/b^2) * diff(diff(p_t,n),n),...
    (2/(a*b)) * diff(diff(p_t,E),n),...
    ];

P = [
    (1/a^2) * diff(diff(p,E),E);...
    (1/b^2) * diff(diff(p,n),n);...
    (2/(a*b)) * diff(diff(p,E),n);...
    ];

% Multiplicacion de [P]^T*[D]*[P] e Integraci?n para obtener matriz de
% rigidez ref. Eq. (8.51b)
pp_k = P_t * D * P;
pp_k_int = int(int(pp_k,E,[-1 1]),n,[-1 1]);


% Matriz evaluaci?n del polinomio ref. Eq. (8.48)
Ae = [
    1   -1   -1   1     1     1   -1    -1     -1   -1    1      1;...
    0    0   1/b  0   -1/b  -2/b   0     1/b   2/b   3/b  -1/b  -3/b;...
    0   -1/a  0   2/a  1/a    0   -3/a  -2/a  -1/a  0     3/a   1/a;...
    
    1    1   -1   1    -1     1    1    -1      1   -1   -1     -1;...
    0    0   1/b  0    1/b  -2/b   0     1/b  -2/b  3/b  1/b    3/b;...
    0   -1/a  0  -2/a  1/a    0   -3/a   2/a  -1/a  0     3/a   1/a;...
    
    1    1    1   1     1     1    1     1      1    1    1      1;...
    0    0   1/b  0    1/b   2/b   0     1/b   2/b   3/b  1/b   3/b;...
    0   -1/a  0  -2/a -1/a    0   -3/a  -2/a  -1/a  0    -3/a  -1/a;...
    
    1   -1    1   1   -1     1   -1     1     -1    1    -1     -1;...
    0    0   1/b  0   -1/b   2/b  0     1/b  -2/b   3/b  -1/b  -3/b;...
    0   -1/a  0   2/a -1/a    0  -3/a   2/a  -1/a   0    -3/a  -1/a;...
    
    ];


% Evaluacion de las anteriores ecuaciones con varibales simbolicas
% usando valores definidos en las variables de netrada. La evaluaci?n
% se hace usando la funci?n 'subs()'

Iz = h^3 / 12;                   % Segundo momento de ?rea rectangular

rho = densidad_rho;
n_nodes = nd;
e = n_nodes - 1;

a = Lx/(2*e);
b = Ly/(2*e);
Em = modulo_young;
v =  coef_poisson;

% Evaluaci?n de Ae y creaci?n de Ae inversa e inversa-transpuesta para
% obtener matrices elementales de inercia y rigidez ref. Eq. (8.51ab)
Ae_subs = subs(Ae);
Ae_inv = inv(Ae_subs);
Ae_inv_t = Ae_inv';

% Calculo de matriz de inercia [Me]
M_e = double(Ae_inv_t * rho * h * pp_m_int * a * b * Ae_inv);

% Calculo de matriz de rigidez [Ke]
pp_k_int_subs = subs(pp_k_int);
K_e = double(Ae_inv_t * Iz * pp_k_int_subs * a * b * Ae_inv);

% Ensamblaje de la matrices globales
Mg = double(assembly_global(M_e,n_nodes));
Kg = double(assembly_global(K_e,n_nodes));


% Eliminar Nodos de los extremos, condiciones de frontera: no hay desplza -
% miento en los bordes de la placa

no_dofs = [];
contador = 1;
for i=1:n_nodes
    no_dofs(contador) = dof_n(i,1,n_nodes);
    no_dofs(contador + 1) = no_dofs(contador) + 1;
    no_dofs(contador + 2) = no_dofs(contador + 1) + 1;
    contador = contador + 3;
    
    no_dofs(contador) = dof_n(i,n_nodes,n_nodes);
    no_dofs(contador + 1) = no_dofs(contador) + 1;
    no_dofs(contador + 2) = no_dofs(contador + 1) + 1;
    contador = contador + 3;
    
    no_dofs(contador) = dof_n(1,i,n_nodes);
    no_dofs(contador + 1) = no_dofs(contador) + 1;
    no_dofs(contador + 2) = no_dofs(contador + 1) + 1;
    contador = contador + 3;
    
    no_dofs(contador) = dof_n(n_nodes,i,n_nodes);
    no_dofs(contador + 1) = no_dofs(contador) + 1;
    no_dofs(contador + 2) = no_dofs(contador + 1) + 1;
    contador = contador + 3;
end

no_dofs = unique(no_dofs);

Kg(:,no_dofs) = [];
Kg(no_dofs,:) = [];

Mg(:,no_dofs) = [];
Mg(no_dofs,:) = [];

% C?lculo del espectro de [ke] y [Me] en un problema de vectores y valores 
% propios en forma general

[A,LAMBDAI] = eig(Kg,Mg);
LAMBDA = LAMBDAI*ones(length(LAMBDAI),1);


% Resultados
disp('*****************************************');
disp(' ');
disp('Modelamiento de los Modos de Flexion')
disp('de una Placa Fina por Elementos Finitos')
disp(' ');
disp('*****************************************');
disp('Caracteristicas de la placa');
disp(strcat('Largo(m): ',num2str(Lx) ))
disp(strcat('Ancho(m): ',num2str(Ly) ))
disp(strcat('grosor(m): ',num2str(h) ))
disp('*****************************************');
disp('Caracteristicas del Material');
disp(strcat('Densidad(kg / m^3): ',num2str(densidad_rho) ))
disp(strcat('Modulo de Young: ',num2str(modulo_young) ))
disp(strcat('Coeficiente de Poisson: ',num2str(coef_poisson) ))
disp('*****************************************');
disp('Detalles del modelamiento');
disp(strcat('Numero total de nodos: ',num2str(n_nodes * n_nodes) ))
disp(strcat('Numero total de elementos: ',num2str(e * e) ))

disp('*****************************************');
disp('Resultados');
% Frecuencias naturales del sistema
Freqs = round(( LAMBDA.^(0.5) )./(2*pi));
% Se imprimen en consola las frecuencias
disp('15 Frecuencias naturales')
disp(Freqs(1:15))



% Grafico de las 12 frecuencias mas bajas

ButtonName = questdlg('Desea graficar los resultados?');

if strcmp(ButtonName, 'Yes')
    
    set(0,'DefaultAxesColor',[1 1 1])
    figure()
    r = 5;
    t = 3;
    for u=1:r*t
        
        n_freq = u;
        C = A(:,n_freq);
        idx = 1:3:length(A);
        C = C(idx);
        
        subplot(r,t,u)
        
        mesh_size = length(C)^0.5;
        [X,Y] = meshgrid(0:1/(mesh_size-1):1,0:1/(mesh_size-1):1);
        surf(X,Y,reshape(C,[mesh_size,mesh_size]))
        view(70,55)
        colormap('jet')
        axis off
        title_str = strcat('',...
            'fn=',num2str(Freqs(n_freq)), ' Hz');
        
        title(title_str)
        set(0,'DefaultAxesColor',[1 1 1])
    end
    
end

ButtonName2 = questdlg('Desea generar una animacion a partir de los resultados(Atencion: El proceso puede tomar mucho tiempo)?');

if strcmp(ButtonName2, 'Yes')
    
    % Crear Animacion
    
    clear T;
    n_frames = 80;
    n_loop = 3;
    counter = 1;
    n_elementos = e^2;
    
    figure;
    for i=1:10
        n_freq = i;
        
        for j= 1:n_loop
            for k = 0:n_frames
                
                
                f = cos((k/n_frames)*(2)*pi);
                C = (f)*A(:,n_freq);
                idx = 1:3:length(A);
                C = C(idx);
                
                mesh_size = length(C)^0.5;
                [X,Y] = meshgrid(0:1/(mesh_size-1):1,0:1/(mesh_size-1):1);
                
                surf(X,Y,reshape(C,[mesh_size,mesh_size]))
                axis([ 0 1 0 1 -5 5]);
                caxis([-5 5])
                view(70 + counter*0.2,32)
                %              view(70,32)
                colormap('jet')
                axis off
                
                cb = colorbar('south');
                title(cb,'Deformacion')
                
                title_str = { strcat('Simulacion con: ',...
                    num2str(e),'x',num2str(e),...
                    ' Elementos');...
                    strcat('Modo #',num2str(n_freq),...
                    '-Frecuencia:',num2str(Freqs(n_freq)), ' Hz')
                    };
                
                dim = [.35 .6 .3 .3];
                an = annotation('textbox',dim,'String',title_str,'FitBoxToText','on',...
                    'BackgroundColor', 'white');
                an.FontSize = 14;
                
                
                T(counter) = getframe(gcf,[40 10 500 400]);
                counter = counter +1;
                
                
                hold off
            end
            
            
        end
    end
    
    vidObj = VideoWriter(strcat('anim_',...
        num2str(n_elementos),'_elementos.mp4'),'MPEG-4');
    open(vidObj);
    writeVideo(vidObj,T);
    close(vidObj);
    
end