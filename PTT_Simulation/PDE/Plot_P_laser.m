%% Ploteo potencia del láser
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