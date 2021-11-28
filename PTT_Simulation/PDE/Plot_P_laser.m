%% Plotting of the laser power:
% Author: Mauricio Cespedes Tenorio
% Date: Nov. 27th, 2021
% Copyright: Laboratorio de Investigacion en Ingenieria Biomedica, UCR. 2021
% Description: I used this code to plot in a prettier way the laser power. Must
% be executed before 'Ploteo_multiple.m'. Two important details are: 
% 1. potencia_L_80: comes from the variable potencia_L from the file 'PTT_controlled.m'. I 
% performed 3 simulations (with T_max = {60,70,80}). Then, I renamed the variable potencia_L
% for each of this results (ex: potencia_L_70 or potencia_L_80). Then I ran this code.
% 2. The variable potencia_laser_80 must be renamed for each T_max (ex: potencia_laser_70 or
% potencia_laser_80).
potencia_laser_80 = [];
indice = 1;
for segundos=0:0.01:699.99
    if(segundos==tiempo(indice))
     %display(indice)
     %display(segundos)
     potencia_actual = potencia_L_80(indice);
     indice = indice + 1;
    end 
    potencia_laser_80 = cat(2,potencia_laser_80, potencia_actual);
end