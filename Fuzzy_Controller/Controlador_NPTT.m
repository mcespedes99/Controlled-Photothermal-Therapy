%% Fuzzy logic controller configuration.
% Author: Mauricio Cespedes Tenorio
% Date: Nov. 28th, 2021
% Copyright: Laboratorio de Investigacion en Ingieria Biomedica, UCR. 2021
% Description: This code implements a Fuzzy Logic Controller with a Mandami Fuzzy Interference
% System. 4 inputs, 1 output and 13 rules are used. If the controller is updated. The file
% ".\MAT_files\fuzzy_controller_NPTT.mat" must be updated before executig 'PTT_controlled.m'.
%

clear;clc;
% Creation of a Mamdami FIS with name "controller"
fis = mamfis('Name',"controller");
%% First input variable: Error Thermal Damage
fis = addInput(fis,[0 100],'Name',"Error Thermal Damage");
% Addition of membership functions:
fis = addMF(fis,"Error Thermal Damage","trimf",[-10 0 10],'Name',"Complete");
fis = addMF(fis,"Error Thermal Damage","trapmf",[0 5 100 120],'Name',"Incomplete");


%% Second input variable: Error Future Maximum Temperature (setpoint max=85°C) Calculado como: e = (setpoint-Temp actual) 
fis = addInput(fis,[-25 50],'Name',"Error Max. Temperature");
% Addition of membership functions:
fis = addMF(fis,"Error Max. Temperature","trapmf",[-30 -25 -5 -2.5],'Name',"Very negative");
fis = addMF(fis,"Error Max. Temperature","trimf",[-5 -2.5 0],'Name',"Negative");
fis = addMF(fis,"Error Max. Temperature","trimf",[-2.5 0 10],'Name',"Zero");
fis = addMF(fis,"Error Max. Temperature","trapmf",[0 10 20 25],'Name',"Positive");
fis = addMF(fis,"Error Max. Temperature","trapmf",[20 25 50 60],'Name',"Very positive");

%% Third input variable: Error Future Maximum Temperature in Healthy Tissue. Setpoint max=50°C. e = (setpoint-Temp actual)
fis = addInput(fis,[-10 15],'Name',"Error Temperature Healthy"); 
% Addition of membership functions:
fis = addMF(fis,"Error Temperature Healthy","trapmf",[-15 -10 -2 -1.5],'Name',"Very negative");
fis = addMF(fis,"Error Temperature Healthy","trapmf",[-2 -1.5 -1 0],'Name',"Negative");
fis = addMF(fis,"Error Temperature Healthy","trapmf",[-1 -0.05 0.05 1.5],'Name',"Zero");
fis = addMF(fis,"Error Temperature Healthy","trapmf",[0 1.5 15 17],'Name',"Positive");

%% Fourth input variable: Rate Thermal Damage in Healthy Tissue
%fis = addInput(fis,[-0.0325 0.075],'Name',"Rate TD Healthy");
% Addition of membership functions:
%fis = addMF(fis,"Rate TD Healthy","trapmf",[-0.0325 0 0.025 0.0325],'Name',"Slow");
%fis = addMF(fis,"Rate TD Healthy","trapmf",[0.025 0.0325 0.075 0.1],'Name',"Fast");

%% Plot of inputs:
figure
plotmf(fis,'input',1);
figure
plotmf(fis,'input',2);
figure
plotmf(fis,'input',3);
%figure
%plotmf(fis,'input',4);

%% Output of Fuzzy Set: delta U en Watts
fis = addOutput(fis,[-0.5 0.5],'Name',"Cambio Potencia Promedio");
fis = addMF(fis,"Cambio Potencia Promedio","trimf",[-0.75 -0.5 -0.25],'Name',"Very negative");
fis = addMF(fis,"Cambio Potencia Promedio","trimf",[-0.5 -0.25 0],'Name',"Negative");
fis = addMF(fis,"Cambio Potencia Promedio","trimf",[-0.25 0 0.25],'Name',"Zero");
fis = addMF(fis,"Cambio Potencia Promedio","trimf",[0 0.25 0.5],'Name',"Positive");
fis = addMF(fis,"Cambio Potencia Promedio","trimf",[0.25 0.5 0.75],'Name',"Very positive");

%Plot of output
figure
plotmf(fis,'output',1)

%% Definition of rules:
ruleList = [2 5 4 5 1 1; % Daño térmico incompleto, e Fut T max muy + y e T sano + => delta U muy +  
            2 4 4 4 1 1; % Daño térmico incompleto, e Fut T max + y e T sano + => delta U +
            2 3 4 3 1 1; % Daño térmico incompleto, e Fut T max 0 y e T sano + => delta U 0
            2 2 4 2 1 1; % Daño térmico incompleto, e Fut T max - y e T sano + => delta U -
            2 5 3 3 1 1; % Daño térmico incompleto, e Fut T max muy + y e T sano 0 => delta U 0
            2 4 3 3 1 1; % Daño térmico incompleto, e Fut T max + y e T sano 0 => delta U 0
            2 3 3 3 1 1; % Daño térmico incompleto, e Fut T max 0 y e T sano 0 => delta U 0
            2 2 3 2 1 1; % Daño térmico incompleto, e Fut T max - y e T sano 0 => delta U -
            2 5 2 2 1 1; % Daño térmico incompleto, e Fut T max muy + y e T sano - => delta U -
            2 4 2 2 1 1; % Daño térmico incompleto, e Fut T max + y e T sano - => delta U -
            2 3 2 2 1 1; % Daño térmico incompleto, e Fut T max 0 y e T sano - => delta U -
            2 2 2 2 1 1; % Daño térmico incompleto, e Fut T max - y e T sano - => delta U -
            1 1 1 1 1 2];% Daño térmico completo, e Fut T max muy - o e T sano muy - => delta U muy -

% Addition of rules to the model
fis = addRule(fis,ruleList);