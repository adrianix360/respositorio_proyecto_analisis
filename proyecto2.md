# Cargar datos:
```matlab
libre = readmatrix('/MATLAB Drive/analisis/vibracion libre.txt');
```

# Cálculo de periodo promedio:
```matlab
t=libre(:,1);
y= libre(:,3);


[pks,locs] = findpeaks(y,t,'MinPeakHeight',4);


primeros_picos = pks(1:3);
primeros_tiempos = locs(1:3);


disp(primeros_tiempos)
```

```matlabTextOutput
    6.8300
    7.3800
    7.9300
```

```matlab


T1 = primeros_tiempos(2)-primeros_tiempos(1);
T2 = primeros_tiempos(3)-primeros_tiempos(2);


T = mean([T1 T2]);
```

# Encontrar amplitud máxima:
```matlab
[ymax,idx] = max(y);


tmax = t(idx);


t_objetivo = tmax + T;


[~,idx2] = min(abs(t-t_objetivo));


y2 = y(idx2);


A1 = ymax;
A2 = y2;
```

# Ploteo de desplazamiento vs tiempo:
```matlab
figure
plot(libre(:,1),libre(:,3),'-b')
legend('Edificio bajo')
xlabel('Tiempo (s)')
ylabel('Desplazamiento (mm)')
title('Respuesta del Edificio Bajo a Vibraciones Libres');
grid on;
```

![figure_0.png](./proyecto2_media/figure_0.png)

```matlab


%Par x-y: 7.38,11.09
%7.93,10.09
%6.83,12.1
%0.55 y 0.55
```

Comienzo calculando la ecuación de la vibración:

x(t)= Xe^{\-\\xi\*\\omega\_n\*t}\*sen(\\omega\_d\*t + \\phi)