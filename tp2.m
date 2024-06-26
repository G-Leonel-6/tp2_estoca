clear all

%tp2 ejercicio 1

% leemos la señal de entrada

x1 = csvread("Archivos/eeg_ojos_abiertos_t7.csv");
N1 = length(x1); %numero de muestras

x2 = csvread("Archivos/eeg_ojos_cerrados_t7.csv");
N2 = length(x2); %numero de muestras

fs = 200;

%a) Defina una funci´on con prototipo ar model(x, P), donde x es la se˜nal y P el orden del
% modelo. La funci´on debe retornar el vector a con los coeficientes y la ganancia G.

% Utilice la funci´on del punto anterior para estimar todos los par´ametros suponiendo ´ordenes
% P = {2, 13, 30}. Grafique en cada caso (para ω ∈ [0, π)) el periodograma superpuesto a la
% PSD del modelo que cumple con la ecuaci´on (6), en base a los par´ametros estimados


%calculo el periodograma de la señal original

SX_1 = (1/N1) * (abs(fft(x1))).^2; %ojos abiertos

SX_2 = (1/N1) * (abs(fft(x2))).^2; %ojos cerrados

%parametros para P = 2
[a2_ojos_abietos, G2_ojos_abiertos] = ar_model(x1, 2); %obtengo los coef G y a_i

%calculo la psd a partir de los parametros obtenidos
[H2_ojos_abiertos, w] = freqz(G2_ojos_abiertos, [1 -a2_ojos_abietos'], 'whole', N1); 
SX2_ojos_abiertos = (abs(H2_ojos_abiertos)).^2; %PSD teorica  

% ojos cerrados
[a2_ojos_cerrados, G2_ojos_cerrados] = ar_model(x2, 2); %obtengo los coef G y a_i

%calculo la psd a partir de los parametros obtenidos
[H2_ojos_cerrados, w] = freqz(G2_ojos_cerrados, [1 -a2_ojos_cerrados'], 'whole', N1); 
SX2_ojos_cerrados = (abs(H2_ojos_cerrados)).^2; %PSD teorica  

%parametros para P = 13
[a13_ojos_abiertos, G13_ojos_abiertos] = ar_model(x1, 13);

%PSD para P = 13
H13_ojos_abiertos = freqz(G13_ojos_abiertos, [1 -a13_ojos_abiertos'], 'whole', N1);
SX13_ojos_abiertos = (abs(H13_ojos_abiertos)).^2; %PSD teorica 

[a13_ojos_cerrados, G13_ojos_cerrados] = ar_model(x2, 13);

%PSD para P = 13
H13_ojos_cerrados = freqz(G13_ojos_cerrados, [1 -a13_ojos_cerrados'], 'whole', N1);
SX13_ojos_cerrados = (abs(H13_ojos_cerrados)).^2; %PSD teorica 

%parametros para P = 30
[a30_abiertos, G30_abiertos] = ar_model(x1, 30);


%PSD para P=30
H30_ojos_abiertos = freqz(G30_abiertos, [1 -a30_abiertos'], 'whole', N1);
SX30_ojos_abiertos = (abs(H30_ojos_abiertos)).^2; %PSD teorica 

[a30_cerrados, G30_cerrados] = ar_model(x2, 30);


%PSD para P=30
H30_ojos_cerrados = freqz(G30_cerrados, [1 -a30_cerrados'], 'whole', N1);
SX30_ojos_cerrados = (abs(H30_ojos_cerrados)).^2; %PSD teorica

f = w*fs;

%grafico los periodogramas
%P=2
figure();
plot(w,10*log10(SX_1), LineWidth = 1); %periodograma de la señal original 
hold on 
plot(w, 10*log10(SX2_ojos_abiertos), "red",LineWidth = 1.5); %PSD para P=2
plot(w, 10*log10(SX13_ojos_abiertos), "green", LineWidth = 1.5); %PSD para P=13
plot(w, 10*log10(SX30_ojos_abiertos), "yellow",LineWidth = 1.5); %PSD para P=30
title("PSD del EEG para ojos abiertos");
legend("Periodograma de la señal original", "PSD teorica con P=2", "PSD teorica con P=13","PSD teorica con P=30",'Location', 'best');
xlabel("$\omega$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
xlim([0 pi]);
grid on;
hold off

% grafico para los datos con ojos cerrados
figure();
plot(w,10*log10(SX_2), LineWidth = 1); %periodograma de la señal original 
hold on 
plot(w, 10*log10(SX2_ojos_cerrados), "red",LineWidth = 1.5); %PSD para P=2
plot(w, 10*log10(SX13_ojos_cerrados), "green", LineWidth = 1.5); %PSD para P=13
plot(w, 10*log10(SX30_ojos_cerrados), "yellow",LineWidth = 1.5); %PSD para P=30
title("PSD del EEG para ojos cerrados");
legend("Periodograma de la señal original", "PSD teorica con P=2", "PSD teorica con P=13","PSD teorica con P=30",'Location', 'best');
xlabel("$\omega$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
xlim([0 pi]);
grid on;
hold off

%============================Ejercicio 2==============================================%
% Aplicando el m´etodo de Welch (Ap´endice 4.1), para una ventana de Hamming de largo
% M = 80 y un solapamiento del 50 %, estime la PSD de ambos registros de EEG. Grafique la
% PSD estimada con Welch junto a la PSD del modelo para P = 13.

% Parámetros del método de Welch
M = 80; % Ancho del segmento
overlap = M / 2; % Solapamiento del 50%


%estimamos la PSD para los datos del ECG con ojos abiertos mediante el
%metodo Welch

Sxx_ojos_abiertos = metodo_Welch(x1, N1, M, overlap);
Sxx_ojos_cerrados = metodo_Welch(x2, N1, M, overlap);
w1 = linspace(0,2*pi,M);
f1 = w1*fs;

% Graficar la PSD para ojos abiertos
figure();
hold on
plot(f,10*log10(SX13_ojos_abiertos), "red", LineWidth = 1); %psd para P=13
plot(f1, 10*log10(Sxx_ojos_abiertos), "blue", LineWidth = 1);
xlabel("$f[Hz]$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
title('PSD usando el método de Welch (ojos abiertos)');
legend("PSD teorica P=13", "PSD estimada", 'Location', 'best');
xlim([0 pi*fs]);
grid on
hold off

% Graficar la PSD para ojos abiertos
figure();
hold on
plot(f,10*log10(SX13_ojos_cerrados), "red", LineWidth = 1); %psd para P=13
plot(f1, 10*log10(Sxx_ojos_cerrados), "blue", LineWidth = 1);
xlabel("$f[Hz]$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
title('PSD usando el método de Welch (ojos cerrados)');
legend("PSD teorica P=13", "PSD estimada", 'Location', 'best');
xlim([0 pi*fs]);
grid on
hold off

%============================EJERCICIO 3======================================
% Utilice los par´ametros obtenidos del modelo para sintetizar las se˜nales EEG artificiales,
% suponiendo un orden P = 13. Luego estime los periodogramas mediante el m´etodo de Welch
% de forma analoga al ejercicio anterior. Compare en un gr´afico la PSD estimada para las se˜nales
% sintetizadas y la PSD del modelo te´orico seg´un la ecuaci´on (6).

%para obtener la señal sintetizada, primero genero ruido gaussiano blanco con N muestras
%luego utilizando filter con los parametros obtenidos por ar_model, obtengo
%la señal sintetizada

ruido = normrnd(0,1,1,N1); %ruido blanco

x13_ojos_abiertos = filter(G13_ojos_abiertos, [1 -a13_ojos_abiertos'], ruido); %señal sintetizada con P=13

%utilizamos la señal sintetizada con P=13 para estimar la PSD mediante el
%metodo Welch

SX13_welch_ojos_abiertos = metodo_Welch(x13_ojos_abiertos', N1, M, overlap);

%ojos cerrados
x13_ojos_cerrados = filter(G13_ojos_cerrados, [1 -a13_ojos_cerrados'], ruido);

SX13_welch_ojos_cerrados = metodo_Welch(x13_ojos_cerrados', N1, M, overlap);

%grafico para ojos abiertos
figure();
hold on
plot(f, 10*log10(SX13_ojos_abiertos), "red", LineWidth = 1); %PSD teorica para P=13
plot(f1, 10*log10(SX13_welch_ojos_abiertos), "blue", LineWidth = 1);
xlabel("$f[Hz]$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
title('PSD señal sintetizada usando el método de Welch (ojos abiertos)');
legend("PSD teorica P=13", "PSD estimada", 'Location', 'best');
xlim([0 pi*fs])
grid on;
hold off

%grafico para ojos cerrados
figure();
hold on
plot(f, 10*log10(SX13_ojos_cerrados), "red", LineWidth = 1); %PSD teorica para P=13
plot(f1, 10*log10(SX13_welch_ojos_cerrados), "blue", LineWidth = 1);
legend("PSD teorica P=13", "PSD estimada", 'Location', 'best');
xlabel("$f[Hz]$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
title('PSD señal sintetizada usando el método de Welch (ojos cerrados)');
xlim([0 pi*fs])
grid on;
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


%estima la PSD de una señal mediante el metodo Welch. Recibe como
%parametros, ademas de la señal, el largo de la señal,
%el numero de segmentos M y el solapamiento
function Pxx = metodo_Welch(x, N, M, overlap)

    % Calcular PSD usando el método de Welch
    K = M - overlap; % Distancia entre segmentos

    % Calculando la PSD usando el método de Welch
    L = floor((N - M) / K) + 1; % Número total de segmentos
    Pxx = zeros(M);

    for i = 1:L
        segment = x((i-1)*K + 1:(i-1)*K + M); % Seleccionar segmento
        %segment = segment - mean(segment); % Remover la media del segmento
        window = hamming(M); % Ventana de Hamming
        segment = segment .* window; % Aplicar ventana
        Pxx = Pxx + abs(fft(segment)).^2; % PSD del segmento
    end

    % Promediando las PSD de los segmentos
    Pxx = (1/N) * Pxx;
    powV = (1/M) * sum((abs(window').^2));
    Pxx = Pxx/powV;

end
