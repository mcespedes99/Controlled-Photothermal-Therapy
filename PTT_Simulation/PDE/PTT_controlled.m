%% Simulation through PDE Toolbox to emulate temperature distribution on tissue under hyperthermia.
% Author: Mauricio Cespedes Tenorio
% Date: Nov. 28th, 2021
% Copyright: Laboratorio de Investigacion en Ingenieria Biomedica, UCR. 2021
% Description: This code implements a simulation through the Partial Differential Equations 
% Toolbox from MATLAB to emulate the temperature distribution in a tissue under laser
% excitation (hyperthermia). To be able to run this code successfully, you must run the
% Montecarlo simulation ('Montecarlo_algorithm.m'), have the controller configuration
% saved under '.\MAT_files\fuzzy_controller_NPTT.mat' and execute the following files first:
% 'Geometry_creation.m' and 'Montecarlo_to_PDE.m'. To plot the results after multiple executions
% refer to the files 'Ploteo_multiple.m' and 'Plot_P_laser.m'. This code is based on the file
% 'PDE_PTT.m'.
%

%% 1. Configuration of the simulation:
clear; clc;
% Simulation time: My recommendation is to execute the file by parts; i.e., not executing it
% from 0s to 700s, but to execute first from 0 to 100s, then from 100s to 200s, and so on.
t_total = 700; %5400s = 1.5 horas
% Intervalos de tiempo a utilizar:
t_iterar = 10; %Cada 5 segundos
P_max = 4; %Potencia máxima del láser

%% 2. Definition of parameters for thermal damage approximation
A_t = 1.98e106;
A_s = 1.18e44;
E_t = 6.67e5;
E_s = 3.02e5;
R_cte = 8.3145;

%% 3. Loading of Fuzzy Logic Controller:
load(".\MAT_files\fuzzy_controller_NPTT.mat");

%% 4. Loading of mat file containing the power density (W/cm^3):
load('.\MAT_files\Heat_rate_laser.mat');


%% 5. Initial conditions: 
% Laser power in Watts:
global flujo
flujo(4,:)=1*flujo(4,:);
% Laser state (ON=1, OFF=0):
global L
L = 0;
%Inital time:
t=t_iterar;

%% 6. Creation of thermal model that is going to be use a the controlled system
% a) Creation of the PDE that we need to solve the problem:
thermalmodel = createpde(1);

% b) Creation of Geometry:
load(".\MAT_files\Geometry_NPTT.mat")
thermalmodel.Geometry = g;

% c) Definition of the mesh:
mesh=generateMesh(thermalmodel,'Hmax',0.25,'Hmin', 0.05);
%figure
%pdemesh(mesh)

% d) Themal properties of tissue:
coef_tejido = specifyCoefficients(thermalmodel,'Cell',1,'m',0,'d',@d,"c",@k,"a",@a_c,"f",@f);
    
% e) Boundary conditions: All internal faces have constant temperature
applyBoundaryCondition(thermalmodel,'dirichlet','Face',2,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',4,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',5,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',6,'u',37);
applyBoundaryCondition(thermalmodel,'dirichlet','Face',1,'u',37);

% Function to find the closest node to an specific coordinate:
getClosestNode = @(p,x,y,z) min((p(1,:) - x).^2 + (p(2,:) - y).^2 + (p(3,:) - z).^2);

%% 7. Variables that are used by the controller:
% Thermal Damage at the coldest point inside the tumor:
TD = [0];
error_TD = [100];
T_TD = [37];
[~,nodo_TD] = getClosestNode(mesh.Nodes, 0.5,0,0);

% Temperature at hottest point in Geometry:
Max_T = [37];
%error_Max_T = [40]; % Not used right now
[~,nodo_Max_T] = getClosestNode(mesh.Nodes, 0,0,0);

%Prediction Max temperature:
error_Max_T_fut = [43];

% Max temperature rate:
% Max_T_rate = [0]; % Not used right now

% Temperature Healthy Tissue:
T_healthy = [37];
[~,nodo_T_H] = getClosestNode(mesh.Nodes, 0,0.7,0);
%error_T_healthy = [11]; % Not used right now

%Prediction Temperature Healthy Tissue:
error_T_healthy_fut = [9];

% Thermal Damage Healthy Tissue:
TD_healthy = [0];


% d(Thermal Damage Healthy Tissue)/dt:
%dt_TD_healthy = [0]; % Not used right now

tiempo = [0];
potencia_L = [0];
%Para controlador Mamdami:
delta_u = 0;

%% Simulación:
while(t<=t_total)
    %Controller action:
    %Controlador v2
    %L = evalfis(fis,[TD(end) Max_T(end) Max_T_rate(end) TD_healthy(end) dt_TD_healthy(end)]); % Not used right now
    %Controlador v3
    %L = evalfis(fis,[TD(end) Max_T(end) Max_T_rate(end) T_healthy(end)]); % Not used right now
    
    %Controlador Mamdami:
    delta_u = evalfis(fis,[error_TD(end) error_Max_T_fut(end) error_T_healthy_fut(end)]);
    if(L+delta_u > P_max)
        L = P_max;
    elseif(L+delta_u < 0)
        L = 0;
    else
        L = L + delta_u;
    end
    coef_tejido.f = @f;
    potencia_L = cat(2, potencia_L, L);
    
    % Initial condition for iteration:
    if(t~=t_iterar) %If time is not 0, the prior initial condition is deleted
        delete(thermalmodel.InitialConditions);
    else
       R = 37; 
    end
    setInitialConditions(thermalmodel,R);
    applyBoundaryCondition(thermalmodel,'dirichlet','Face',2,'u',37);
    applyBoundaryCondition(thermalmodel,'dirichlet','Face',4,'u',37);
    applyBoundaryCondition(thermalmodel,'dirichlet','Face',5,'u',37);
    applyBoundaryCondition(thermalmodel,'dirichlet','Face',6,'u',37);
    applyBoundaryCondition(thermalmodel,'dirichlet','Face',1,'u',37);
    
    % Time used for the PDE solution: 
    tlist = [0:5:t_iterar];
    
    %Solution of the PDE:
    fprintf("\n Tiempo de simulación: ");
    disp(t);
    R = solvepde(thermalmodel,tlist);
    T = R.NodalSolution;
    
    % Temperature at point used for the Thermal Damage Approximation:
    T_TD = cat(2, T_TD, T(nodo_TD,end));
    Temp_TD = T(nodo_TD,:);
    if(T(nodo_TD,1)>=40)
      TD_actual = 100*(1-exp(-trapz(tlist,A_t*exp(-1*E_t./(R_cte*(Temp_TD+273.15))))));
    else
      TD_actual = 0;    
    end
    TD_actual = TD_actual + TD(end);
    if(TD_actual>=100)
        TD_actual = 100;
    end
    TD = cat(2, TD, TD_actual);
    error_TD = cat(2, error_TD, 100-TD(end));
    
    % Temperature at hottest point in Geometry:
    actual_Max_T = T(nodo_Max_T,end);
    Max_T = cat(2, Max_T, actual_Max_T);
    %error_Max_T = cat(2, error_Max_T, 75-Max_T(end)); % Not used right now
    delta_Max_T = Max_T(end)-Max_T(end-1);
    Max_T_fut = Max_T(end)+delta_Max_T;
    % Actualizar error de T futuro:
    error_Max_T_fut = cat(2, error_Max_T_fut, 80-Max_T_fut);
    
    % Max Temp rate: % Not used right now
    %dt_Max_T = (Max_T(end)-Max_T(end-1))/t_iterar;
    %Max_T_rate = cat(2, Max_T_rate, dt_Max_T);
    
    % Temperature Healthy Tissue:
    actual_T_H = T(nodo_T_H,end);
    T_healthy = cat(2, T_healthy, actual_T_H);
    %error_T_healthy = cat(2, error_T_healthy, 46-T_healthy(end)); % Not used right now
    delta_T_h = T_healthy(end) - T_healthy(end-1);
    T_healthy_fut = T_healthy(end)+delta_T_h;
    error_T_healthy_fut  = cat(2, error_T_healthy_fut, 46-T_healthy_fut);
    
    % Thermal Damage Healthy Tissue:
    Temp_TDHT = T(nodo_T_H,:);
    if(T(nodo_T_H,1)>=40)
      TDHT_actual = 100*(1-exp(-trapz(tlist,A_s*exp(-1*E_s./(R_cte*(Temp_TDHT+273.15))))));
    else
      TDHT_actual = 0;    
    end
    TDHT_actual = TDHT_actual + TD_healthy(end);
    if(TDHT_actual>=100)
        TDHT_actual = 100;
    end
    TD_healthy = cat(2, TD_healthy, TDHT_actual);
    
    
    % d(Thermal Damage Healthy Tissue)/dt: % Not used right now
    %dt_TDHT = (TD_healthy(end)-TD_healthy(end-1))/t_iterar;
    %t_TD_healthy = cat(2, dt_TD_healthy, dt_TDHT);
    
    fprintf("\n Potencia láser:");
    disp(L);
    fprintf("\n Temperatura TD:");
    disp(Temp_TD(end));
    fprintf("\n Thermal Damage:");
    disp(TD_actual);
    fprintf("\n Thermal Damage Error:");
    disp(error_TD(end));
    fprintf("\n Temperatura Max:");
    disp(actual_Max_T);
    fprintf("\n Temperatura Max Error Fut:");
    disp(error_Max_T_fut(end));
    %fprintf("\n Rate Temperatura Max:"); % Not used right now
    %disp(dt_Max_T); % Not used right now
    fprintf("\n Temperatura TH:");
    disp(actual_T_H);
    fprintf("\n TD HT:");
    disp(TDHT_actual);
    fprintf("\n Error Fut Temp HT:");
    disp(error_T_healthy_fut(end));
    %fprintf("\n d(TD)/dt HT:"); % Not used right now
    %disp(dt_TDHT); % Not used right now
    
    % Update of time counter:
    tiempo = cat(2, tiempo, t);
    t = t+t_iterar;
end
%save(NPTT_con_control);
%Cálculo del daño térmico en tejido sano:
%TD_sano = [0];
%for i=2:81
%    if(T_healthy(i-1)>=40)
%        valor_temp = [T_healthy(i-1) T_healthy(i)];
%        TD_sano(i) = TD_sano(end) + 100*(1-exp(-trapz([0 t_iterar],(1.18e44)*exp(-1*(3.02e5)./(R_cte*(valor_temp+273.15))))));
%    else
%        TD_sano(i)=0;
%    end
%end

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
%display(state.time);
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
