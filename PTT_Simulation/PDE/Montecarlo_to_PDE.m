%% ---------------------------Acople con Montecarlo-------------------------
%Limpieza de consola y del workspace:
clear; clc;

%Carga de matriz de interés:
load('.\MAT_files\Fluence_result_Montecarlo.mat');
%En cada fila de H están las coordenadas de cada tetrahedro. ValoMC usa vectores para definir el tetraedro entonces cada fila de 
%H contiene las filas del array "r" en las que están contenidas los vectores que forman el tetrahedro. En r, lo que están son las
%coordenadas de estos vectores que forman los tetraedros. https://inverselight.github.io/ValoMC/structures3d.html
H = vmcmesh.H;
r = vmcmesh.r;

%Definición del tamaño de matriz de interés: El número de filas de H es
%equivalente al número de tetrahedros o subdivisiones de la geometria.
[m,n] = size(H);

%Variable flujo para guardar el flujo y sus coordenadas. Tiene 4 filas
%porque, para cada elemento de la geometría se van a guardar sus
%coordenadas (x,y,z) y la tasa de calor Q_L asociado a ese elemento.
global flujo;
flujo = zeros(4,m);

%Bucle para recorrer todas las filas y columnas de H:
%El for está pensado para que vaya fila por fila de H (tetraedro por tetraedro),
%sacando las coordenadas de cada uno de los vectores que forman el
%tetraedro. Luego se hace un promedio en cada uno de los ejes x,y,z para
%lograr calcular una coordenada promedio del centro del tetraedro. Luego
%estas coordenadas las pase a esfericas para calcular la tasa de calor
%generada por el laser, dada por la irradiancia (fluence) por mu_a.
%Finalmente, las guarde en el array flujo que creé antes.
for i=1:m
    coordenadas = zeros(4,3);
    contador = 1;
    for j=1:n
       fila_r = H(i,j);
       coordenadas(contador,:)= r(fila_r,:);
       contador = contador + 1;
    end
    xc = mean(coordenadas(:,1))/10;
    yc = mean(coordenadas(:,2))/10;
    zc = mean(coordenadas(:,3))/10;
    [~,~,radio] = cart2sph(xc,yc,zc);
    if(radio <= 0.5)
        flujo_Ne = solution.element_fluence(i)*0.0895438*1000;
    elseif(xc<=0.2)
        flujo_Ne = solution.element_fluence(i)*0.0207638*1000;
    elseif(xc <= 0.4) && (xc > 0.2)
         flujo_Ne = solution.element_fluence(i)*0.008352283*1000;
    elseif(xc > 0.4) && (radio > 0.5)
         flujo_Ne = solution.element_fluence(i)*0.021433576*1000;  
    end
    flujo(:,i)=[xc,yc,zc,flujo_Ne];
end
flujo(1,:) =  round(flujo(1,:),4);
flujo(2,:) =  round(flujo(2,:),4);
flujo(3,:) =  round(flujo(3,:),4);
flujo(4,:) =  round(flujo(4,:),4);
for i = 1:m
    if(flujo(4,i) < 0.001)
       flujo(4,i) = 0;
    end
end

%Esto es para plotearlo bonito. Si quiere entenderlo, me escribe y le
%ayudo, si no logra sacar.
[x,y]=meshgrid(0:0.025:3, -3:0.025:3);
[a,b]=size(x);
Q_XY = zeros(size(x));
Z = 0;
for i=1:b
    actual_x = x(1,i);
    for j=1:a
        actual_y = y(j,1);
        [~,nodo]= min((flujo(1,:) - actual_x).^2 + (flujo(2,:) - actual_y).^2 + (flujo(3,:) - Z).^2);
        Q_XY(j,i)= flujo(4,nodo);
    end
end
