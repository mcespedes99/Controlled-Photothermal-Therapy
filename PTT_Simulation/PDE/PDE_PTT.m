% Ejemplo matlab esfera dentro de cubo: (0.5,0,0), (0,0,0), (0,0.7,0)
%% ---------------------------Montecarlo results for light distribution inside the tissue-------------------------
clear; clc;

%Loading of mat file containing the power density (W/cm^3):
load('tasa_calor_superficial.mat');
%Definir potencia del laser: (actual = 2 W)
global flujo
flujo(4,:)=4*flujo(4,:);
%% -------------------------------PDE Solution------------------------------
% 1) Creation of the PDE's that we need to solve the problem:
%Main PDE
thermalmodel = createpde(1);

%% 2) Creation of Geometry:
load("Geometry_NPTT.mat")
thermalmodel.Geometry = g;

%% 3) Plot of Cells y Faces of the Geometry
%figure('Position',[10,10,800,400]);
%subplot(1,2,1)
%pdegplot(thermalmodel,'FaceAlpha',0.25,'CellLabel','on')
%title('Geometry with Cell Labels')
%subplot(1,2,2)
%pdegplot(thermalmodel,'FaceAlpha',0.25,'FaceLabel','on')
%title('Geometry with Face Labels')

%% 4) Definition of the mesh:
mesh=generateMesh(thermalmodel,'Hmax',0.25,'Hmin', 0.05);
%figure
%pdemesh(mesh)

%% 5) Themal properties needed for the simulation to run:
%Laser on (L=1) or Laser off (L=0)
global L
L = 1;
%Healthy tissue:
specifyCoefficients(thermalmodel,'Cell',1,'m',0,'d',@d,"c",@k,"a",@a_c,"f",@f);
%Damage tissue:
%d_tumor_cte = 0.001*3500;
%k_tumor_cte = 0.00642;
%specifyCoefficients(thermalmodel,'Cell',2,'m',0,'d',d_tumor_cte,"c",k_tumor_cte,"a",@a_tumor,"f",@f_tumor);

%% 6) Boundary conditions Laser:
applyBoundaryCondition(thermalmodel,'dirichlet','Face',2,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',4,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',5,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',6,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',1,'u',37);
%applyBoundaryCondition(thermalmodel,'dirichlet','Face',1,'h',@h,'r',@r);

%Function to find the closest node to an specific coordinate:
getClosestNode = @(p,x,y,z) min((p(1,:) - x).^2 + (p(2,:) - y).^2 + (p(3,:) - z).^2);

% iii) Initial conditions:
%Initial temperature:
R = 37;
setInitialConditions(thermalmodel,R);


%% 7) PDE solution:
%thermalmodel.SolverOptions.MinStep = 0.5;
%model.SolverOptions.AbsoluteTolerance = 1.0e-4;
tlist = [0:5:1000];
fprintf("solución");
R = solvepde(thermalmodel,tlist);
T = R.NodalSolution;

figure;
pdeplot3D(thermalmodel,'ColorMapData',T(:,end)) %Índice 2 = 10 segundos
title(['Temperature at Time 1000 s']);

%% 8) Graficar aumento temporal temperatura en un punto (x,y,z) determinado:
figure
[~,nid2] = getClosestNode(mesh.Nodes, 0,0,-0.5);
plot(tlist, T(nid2,:));
hold on

[~,nid3] = getClosestNode(mesh.Nodes, 0,0,0.5);
plot(tlist, T(nid3,:));

[~,nid4] = getClosestNode(mesh.Nodes, 0,-0.5,0);
plot(tlist, T(nid4,:));

[~,nid5] = getClosestNode(mesh.Nodes, 0,0.5,0);
plot(tlist, T(nid5,:));

[~,nid6] = getClosestNode(mesh.Nodes, 0.5,0,0);
plot(tlist, T(nid6,:));
legend({'T at (0,0,-0.5)','T at (0,0,0.5)','T at (0,-0.5,0)','T at (0,-0.5,0)','T at (0.5,0,0)'});
title 'Temperature of Sphere as a Function of Time';
xlabel 'Time, seconds'
ylabel 'Temperature, degrees-Celsius'
hold off

%Llamado de función para encontrar el nodo más cercano al punto (0,0,0)
%[~,nid1] = getClosestNode(mesh.Nodes, 0,0,0);
%figure;
%plot(tlist, T(nid1,:));
%hold on
%[~,nid7] = getClosestNode(mesh.Nodes, -1,0,0);
%plot(tlist, T(nid7,:));

%[~,nid8] = getClosestNode(mesh.Nodes, -1.1,0,0);
%plot(tlist, T(nid8,:));

%[~,nid9] = getClosestNode(mesh.Nodes, -1.2,0,0);
%plot(tlist, T(nid9,:));

%[~,nid10] = getClosestNode(mesh.Nodes, -1.3,0,0);
%plot(tlist, T(nid10,:));
%grid on
%legend({'T at (0,0,0)', 'T at (-1,0,0)','T at (-1.1,0,0)','T at (-1.2,0,0)','T at (-1.3,0,0)'});
%title 'Temperature of Sphere as a Function of Time';
%xlabel 'Time, seconds'
%ylabel 'Temperature, degrees-Celsius'
%hold off

figure
[~,nid11] = getClosestNode(mesh.Nodes, 0.4,0.3,0);
plot(tlist, T(nid11,:));
hold on

[~,nid12] = getClosestNode(mesh.Nodes, 0.4,-0.3,0);
plot(tlist, T(nid12,:));

[~,nid13] = getClosestNode(mesh.Nodes, 0.4,0,0.3);
plot(tlist, T(nid13,:));

[~,nid14] = getClosestNode(mesh.Nodes, 0.4,0,-0.3);
plot(tlist, T(nid14,:));

[~,nid15] = getClosestNode(mesh.Nodes, 1,0,0);
plot(tlist, T(nid15,:));
legend({'T at (0.4,0.3,0)','T at (0.4,-0.3,0)','T at (0.4,0,0.3)','T at (0.4,0,-0.3)','T at (1,0,0)'});
title 'Temperature of Sphere as a Function of Time';
xlabel 'Time, seconds'
ylabel 'Temperature, degrees-Celsius'
hold off

[x,y]=meshgrid(0:0.05:3, -3:0.05:3);
[a,b]=size(x);
T_XY = zeros(size(x));
Z = 0;
for i=1:b
    actual_x = x(1,i);
    for j=1:a
        actual_y = y(j,1);
        [~,nodo]= getClosestNode(mesh.Nodes, actual_x,actual_y,Z);
        T_XY(j,i)= T(nodo,end);
    end
end
figure
hold on
surf(x,y,T_XY);
xlabel('Profundidad x [cm]');
ylabel('Longitud y [cm]');
zlabel('Temperatura [°C]');

cc = colorbar;                       
c.Label.String = 'Temperatura [(°C)]';


%% Tissue coefficients
%1. Constante "d": Cp*rho
function d = d(location,~)
d = zeros(size(location.x));    
%Para tumor:
[~,~,radio]=cart2sph(location.x,location.y,location.z);
%Para el tumor:
ids_tumor = radio <= 0.5;
d(ids_tumor) = 0.001*3500;
%Para la piel:
ids_piel = location.x <= 0.2 & radio > 0.5;
d(ids_piel) = 0.0012*3800;
%Para adiposo:
ids_adipose = location.x <= 0.4 & location.x > 0.2 & radio > 0.5;
d(ids_adipose) = 0.0023*850;
%Para músculo:
ids_muscle = location.x > 0.4 & radio > 0.5;
d(ids_muscle) = 0.00127*3500;
end

%2. Conductividad térmica: constante "c" de tejido sano
function K = k(location,~)
K = zeros(size(location.x));    
%Para tumor:
[~,~,radio]=cart2sph(location.x,location.y,location.z);
%Para el tumor:
ids_tumor = radio <= 0.5;
K(ids_tumor) = 0.00642;

%Para la piel:
ids_piel = location.x <= 0.2 & radio > 0.5;
K(ids_piel) = 0.0053;

%Para adiposo:
ids_adipose = location.x <= 0.4 & location.x > 0.2 & radio > 0.5;
K(ids_adipose) = 0.0016;

%Para músculo:
ids_muscle = location.x > 0.4 & radio > 0.5;
K(ids_muscle) = 0.0053;
end

%3. Constante "a" tejido sano: -Cp*W(T)
function coef_a = a_c(location,state)
coef_a = zeros(size(location.x));    
%Para tumor:
[~,~,radio]=cart2sph(location.x,location.y,location.z);
%Para el tumor:
ids_t_menor_37 = state.u<37 & radio <= 0.5;
coef_a(ids_t_menor_37) = (3500*(0.833))/1000000;
ids_t_entre_37_42 = state.u>=37 & state.u<=42 & radio <= 0.5;
coef_a(ids_t_entre_37_42) = (3500*(0.833*(state.u(ids_t_entre_37_42)-37).^4.8 / 5438))/1000000;
ids_t_mayor_42 = state.u>42 & radio <= 0.5;
coef_a(ids_t_mayor_42) = (3500*0.416)/1000000;

%Para la piel:
ids_piel_menor_44 = location.x <= 0.2 & radio > 0.5 & state.u<=44;
coef_a(ids_piel_menor_44) = (3500*0.45*(1+9.2*exp(-(state.u(ids_piel_menor_44)-44).^2 / 10)))/1000000;
ids_piel_mayor_44 = location.x <= 0.2 & radio > 0.5 & state.u>44;
coef_a(ids_piel_mayor_44) = (3500*0.45*10.2)/1000000;

%Para adiposo:
ids_adipose_menor_45 = location.x <= 0.4 & location.x > 0.2 & radio > 0.5 & state.u<=45;
coef_a(ids_adipose_menor_45) = (3500*(0.36+0.36*exp(-(state.u(ids_adipose_menor_45)-45).^2 / 12)))/1000000;
ids_adipose_mayor_45 = location.x <= 0.4 & location.x > 0.2 & radio > 0.5 & state.u>45;
coef_a(ids_adipose_mayor_45) = (3500*0.72)/1000000;

%Para músculo:
ids_muscle_menor_45 = location.x > 0.4 & radio > 0.5 & state.u<=45;
coef_a(ids_muscle_menor_45) = (3500*(0.45+3.55*exp(-(state.u(ids_muscle_menor_45)-45).^2 / 12)))/1000000;
ids_muscle_mayor_45 = location.x > 0.4 & radio > 0.5 & state.u>45;
coef_a(ids_muscle_mayor_45) = (3500*4)/1000000;
end

%3. Coeficiente "f": QL+QM-Cp*Tb*W(T)
function f_coef = f(location,state)
global flujo;
global L;
f_coef = nan(size(location.x));
[~,c_s] = size(location.x);
x_array_s = [location.x];
y_array_s = [location.y];
z_array_s = [location.z];
t_array_s = [state.u];
display(state.time);
%cond = isnan(state.time);
%if(cond == 1)
%    f_coef_s = nan(size(location.x));
%else
for i=1:c_s
    x = x_array_s(i);
    y = y_array_s(i);
    z = z_array_s(i);
    t = t_array_s(i);
    [~,~,radio]=cart2sph(x,y,z);
    
    [~,nodo]= min((flujo(1,:) - x).^2 + (flujo(2,:) - y).^2 + (flujo(3,:) - z).^2); 
    P_laser = L*flujo(4,nodo);
    %Primero, se coloca la convección con el ambiente en zona expuesta del tumor:
    if(x == 0) && (radio <= 0.5)
        % Convección Q_c = -h*(Tamb-T) = -5 W/(K*m²)*(25-T) = 
        % modelo_ecuaciones.pdf= Thermal dosage investigation for optimal
        % temperature...-Yatao Ren
        if(t<37)
            f_coef(i)= P_laser+0.001091+(3500*37*(0.833))/1000000 - 0.0005*(25-t);
        elseif(t>=37) && (t<=42)
            f_coef(i)= P_laser+0.001091+(3500*37*(0.833*(t-37)^4.8 / 5438))/1000000 - 0.0005*(25-t);
        else
            f_coef(i)= P_laser+0.001091+(3500*(37)*0.416)/1000000 - 0.0005*(25-t);
        end
    
    %Piel expuesta al ambiente:
    elseif(x == 0) && (radio > 0.5)
        if(t <= 44)
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*0.45*(1+9.2*exp(-(t-44)^2)))/1000000 - 0.0005*(25-t);
        else
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*0.45*10.2)/1000000 - 0.0005*(25-t);
        end
        
    %Para el tumor no expuesto:
    elseif(radio <= 0.5)
        if(t<37)
            f_coef(i)= P_laser+0.001091+(3500*37*(0.833))/1000000;
        elseif(t>=37) && (t<=42)
            f_coef(i)= P_laser+0.001091+(3500*37*(0.833*(t-37)^4.8 / 5438))/1000000;
        else
            f_coef(i)= P_laser+0.001091+(3500*(37)*0.416)/1000000;
        end
    
    %Para la piel:
    elseif(x <= 0.2)
        if(t <= 44)
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*0.45*(1+9.2*exp(-(t-44)^2)))/1000000;
        else
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*0.45*10.2)/1000000;
        end
    
    %Para el tejido adiposo:    
    elseif(x <= 0.4) && (x > 0.2)
        if(t <= 45)
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*(0.36+0.36*exp(-(t-45)^2/12)))/1000000;
        else
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*0.72)/1000000;
        end
    
    %Para el músculo:    
    else
        if(t <= 45)
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*(0.45+3.55*exp(-(t-45)^2 / 12)))/1000000;
        else
           f_coef(i) = P_laser + 0.001091 + (3500*(37)*4)/1000000;
        end
    end
end
end
