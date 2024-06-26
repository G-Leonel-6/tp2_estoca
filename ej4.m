rp = 0.01;
rs = 0.001;
M = 10000;
fs = 200;

x1 = csvread("Archivos/eeg_ojos_abiertos_t7.csv");
x2 = csvread("Archivos/eeg_ojos_cerrados_t7.csv");
N_x1 = length(x1);
N_x2 = length(x2);

N = 280; % 255
F = [0 3*2/fs 5*2/fs 1];
A = [1 1 0 0];
V = [1/rp, 1/rs];

N1 = 280; %128
F1 = [0 3*2/fs 5*2/fs 8*2/fs 10*2/fs 1];
A1 = [0 0 1 1 0 0];
V1 = [1/rs 1/rp 1/rs];

N2 = 280; %128
F2 = [0 8*2/fs 10*2/fs 13*2/fs 15*2/fs 1];
A2 = [0 0 1 1 0 0];
V2 = [1/rs 1/rp 1/rs];

N3 = 270; %128
F3 = [0 13*2/fs 15*2/fs 29*2/fs 31*2/fs 1];
A3 = [0 0 1 1 0 0];
V3 = [1/rs 1/rp 1/rs];

N4 = 280; %255
F4 = [0 29*2/fs 31*2/fs 1];
A4 = [0 0 1 1];
V4 = [1/rs 1/rp];

h = firpm(N, F, A, V);
H = fft(h, M);

h1 = firpm(N1, F1, A1, V1);
H1 = fft(h1, M);

h2 = firpm(N2, F2, A2, V2);
H2 = fft(h2, M);

h3 = firpm(N3, F3, A3, V2);
H3 = fft(h3, M);

h4 = firpm(N4, F4, A4, V4);
H4 = fft(h4, M);

w = linspace(0, fs, M);

% Graficos de filtros

figure
hold on
grid on
title('Filtro de ondas Delta')
xlim([0, 45])
plot(w, abs(H), color = 'blue', LineWidth = 1)
line([0 3], [1.01 1.01], color = 'red', linestyle = '--', LineWidth = 1)
line([0 3], [0.99 0.99], color = 'red', linestyle = '--')
line([5 10], [0.001 0.001], color = 'red', linestyle = '--')
line([5 10], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{D}(f)|')
legend("Respuesta en frecuencia del filtro", "Restricción de ripple", 'location', 'best')
hold off

figure
hold on
grid on
title('Filtro de ondas Theta')
xlim([0, 45])
plot(w, abs(H1), color = 'blue', LineWidth = 1)
line([5 8], [1.01 1.01], color = 'red', linestyle = '--', LineWidth = 1)
line([5 8], [0.99 0.99], color = 'red', linestyle = '--')
line([0 3], [0.001 0.001], color = 'red', linestyle = '--')
line([0 3], [-0.001 -0.001], color = 'red', linestyle = '--')
line([10 12], [0.001 0.001], color = 'red', linestyle = '--')
line([10 12], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{T}(f)|')
legend("Respuesta en frecuencia del filtro", "Restricción de ripple", 'location', 'best')

hold off

figure
hold on
grid on
title('Filtro de ondas Alpha')
xlim([0, 45])
plot(w, abs(H2), color = 'blue', LineWidth = 1)
line([10 13], [1.01 1.01], color = 'red', linestyle = '--', LineWidth = 1)
line([10 13], [0.99 0.99], color = 'red', linestyle = '--')
line([0 8], [0.001 0.001], color = 'red', linestyle = '--')
line([0 8], [-0.001 -0.001], color = 'red', linestyle = '--')
line([15 20], [0.001 0.001], color = 'red', linestyle = '--')
line([15 20], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{A}(f)|')
legend("Respuesta en frecuencia del filtro", "Restricción de ripple", 'location', 'best')
hold off

figure
hold on
grid on
title('Filtro de ondas Beta')
xlim([0, 45])
plot(w, abs(H3), color = 'blue', LineWidth = 1)
line([15 29], [1.01 1.01], color = 'red', linestyle = '--', LineWidth = 1)
line([15 29], [0.99 0.99], color = 'red', linestyle = '--')
line([0 13], [0.001 0.001], color = 'red', linestyle = '--')
line([0 13], [-0.001 -0.001], color = 'red', linestyle = '--')
line([31 35], [0.001 0.001], color = 'red', linestyle = '--')
line([31 35], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{B}(f)|')

legend("Respuesta en frecuencia del filtro", "Restricción de ripple", 'location', 'best')
hold off

figure
hold on
grid on
title('Filtro de ondas Gamma')
xlim([0, 45])
plot(w, abs(H4), color = 'blue', LineWidth = 1)
line([31 45], [1.01 1.01], color = 'red', linestyle = '--', LineWidth = 1)
line([31 45], [0.99 0.99], color = 'red', linestyle = '--')
line([0 29], [0.001 0.001], color = 'red', linestyle = '--')
line([0 29], [-0.001 -0.001], color = 'red', linestyle = '--')
legend("Respuesta en frecuencia del filtro", "Restricción de ripple", 'location', 'best')
xlabel('f[Hz]')
ylabel('|H_{G}(f)|')
hold off

%------------------------------------------***************---------------------------------------------

%sintesis de señales

ruido = normrnd(0,1,1,N); 

[a_13_1, G_13_1] = ar_model(x1, 13); %ojos abiertos

u_13_1 = ruido;
x_13_1 = filter(G_13_1, [1 -a_13_1'], u_13_1);

[a_13_2, G_13_2] = ar_model(x2, 13);   %ojos cerrados

u_13_2 =  ruido;
x_13_2 = filter(G_13_2, [1 -a_13_2'], u_13_2);

%filtrado de señales reales

y_1 = conv(h,x1);
y1_1 = conv(h1, x1);
y2_1 = conv(h2, x1);
y3_1 = conv(h3, x1);
y4_1 = conv(h4, x1);

y_2 = conv(h,x2);
y1_2 = conv(h1, x2);
y2_2 = conv(h2, x2);
y3_2 = conv(h3, x2);
y4_2 = conv(h4, x2);

%filtrado de las señales sintetizadas

y_13_1 = conv(h,x_13_1);
y1_13_1 = conv(h1, x_13_1);
y2_13_1 = conv(h2, x_13_1);
y3_13_1 = conv(h3, x_13_1);
y4_13_1 = conv(h4, x_13_1);

y_13_2 = conv(h,x_13_2);
y1_13_2 = conv(h1, x_13_2);
y2_13_2 = conv(h2, x_13_2);
y3_13_2 = conv(h3, x_13_2);
y4_13_2 = conv(h4, x_13_2);

%calculo de varianza delas señales reales

var_y_1 = [var(y_1) var(y1_1) var(y2_1) var(y3_1) var(y4_1)];
var_y_2 = [var(y_2) var(y1_2) var(y2_2) var(y3_2) var(y4_2)];

%calculo de varianza de señales sintetizadas

var_y_13_1 = [var(y_13_1) var(y1_13_1) var(y2_13_1) var(y3_13_1) var(y4_13_1)];
var_y_13_2 = [var(y_13_2) var(y1_13_2) var(y2_13_2) var(y3_13_2) var(y4_13_2)];

% grafico de potencias

figure
hold on
grid on
xticks([1 2 3 4 5])
plot(var_y_1, "blue", LineWidth = 1);
plot(var_y_13_1, "red", LineWidth = 1);
xticklabels({'D', 'T', 'A', 'B', 'G'});
legend("Potencia de la señal real", "Potencia de la señal sintetizada")
title('Potencia de las señales de ojos abiertos')
hold off


figure
hold on
grid on 
title('Potencia de las señales de ojos cerrados')
xticks([1 2 3 4 5])
plot(var_y_2, 'blue', LineWidth = 1);
plot(var_y_13_2, 'red', LineWidth = 1);
legend("Potencia de la señal real", "Potencia de la señal sintetizada")
xticklabels({'D', 'T', 'A', 'B', 'G'});

hold off

function [a, G] = ar_model(x, P)

    N = length(x); %largo del vector x

    %calculamos la funcion de autocorrelacion de x
    R_x = xcorr(x, 'biased');

    %calculamos la matriz R

    R = zeros(P, P);

    for i = 0:P-1

        R(i+1, :) = R_x(N-i:N+P-1-i);

    end

    %defino el vector r
    r = R_x(N+1:N+P);

    %estimo los coeficientes de a
    a = R^(-1) * r;

    %estimamos la ganancia

    G = (R_x(N) - sum( a'*R_x(N+1:N+P) ))^(1/2);

end
