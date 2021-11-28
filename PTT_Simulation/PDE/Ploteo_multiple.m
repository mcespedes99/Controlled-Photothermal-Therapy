%% Plotting of variables after multiple executions of the file 'PTT_controlled.m'.
% Author: Mauricio Cespedes Tenorio
% Date: Nov. 28th, 2021
% Copyright: Laboratorio de Investigacion en Ingieria Biomedica, UCR. 2021
% Description: As described in the README file for this repository. I performed three simulations
% with different values of T_max. Then, to plot each result, I renamed the variables extracted
% from the execution of the file 'PTT_controlled.m' and plot them. Important: the file
% 'Ploteo_P_laser.m' must be executed before this one.
%

figure
hold on
plot(0:0.01:699.99,potencia_laser_60,'k');
plot(0:0.01:699.99,potencia_laser_70,'k--');
plot(0:0.01:699.99,potencia_laser_80,'k:');
legend("Setpoint T max = 60°C","Setpoint T max = 70°C","Setpoint T max = 80°C")
ylabel("Laser power (W)");
xlabel("Time (s)");
hold off

figure
hold on
plot(tiempo,TD_60,'k');
plot(tiempo,TD_70,'k--');
plot(tiempo,TD_80,'k:');
legend("Setpoint T max = 60°C","Setpoint T max = 70°C","Setpoint T max = 80°C")
ylabel("Tumor thermal damage (%)");
xlabel("Time (s)");
hold off

%figure
%hold on
%plot(tiempo,T_TD_60,'k');
%plot(tiempo,T_TD_70,'k--');
%plot(tiempo,T_TD_80,'k:');
%legend("Setpoint T max = 60°C","Setpoint T max = 70°C","Setpoint T max = 80°C")
%ylabel("Mínima temperatura dentro del tumor (°C)");
%xlabel("Tiempo (s)");
%hold off

figure
hold on
plot(tiempo,Max_T_60,'k');
plot(tiempo,Max_T_70,'k--');
plot(tiempo,Max_T_80,'k:');
legend("Setpoint T max = 60°C","Setpoint T max = 70°C","Setpoint T max = 80°C")
ylabel("Tumor highest temperature (°C)");
xlabel("Time (s)");
hold off

%figure
%plot(tiempo,Max_T_rate);
%ylabel("Pendiente máxima temperatura dentro del tumor (°C/s)");
%xlabel("Tiempo (s)");

figure
hold on
plot(tiempo,T_healthy_60,'k');
plot(tiempo,T_healthy_70,'k--');
plot(tiempo,T_healthy_80,'k:');
legend("Setpoint T max = 60°C","Setpoint T max = 70°C","Setpoint T max = 80°C")
ylabel("Highest healthy tissue temperature (°C)");
xlabel("Time (s)");
hold off

figure
hold on
plot(tiempo,TD_healthy_60,'k');
plot(tiempo,TD_healthy_70,'k--');
plot(tiempo,TD_healthy_80,'k:');
legend("Setpoint T max = 60°C","Setpoint T max = 70°C","Setpoint T max = 80°C")
ylabel("Healthy tissue thermal damage (%)");
xlabel("Tiempo (s)");
hold off

%figure
%plot(tiempo,dt_TD_healthy);
%ylabel("Pendiente daño térmico en el tejido sano circundante (%/s)");
%xlabel("Tiempo (s)");