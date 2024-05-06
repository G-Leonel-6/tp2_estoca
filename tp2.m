clear all

%tp2 ejercicio 1

% leemos la señal de entrada

x = csvread("Archivos/eeg_ojos_abiertos_t7.csv");
N = length(x); %numero de muestras


%a) Defina una funci´on con prototipo ar model(x, P), donde x es la se˜nal y P el orden del
% modelo. La funci´on debe retornar el vector a con los coeficientes y la ganancia G.

% Utilice la funci´on del punto anterior para estimar todos los par´ametros suponiendo ´ordenes
% P = {2, 13, 30}. Grafique en cada caso (para ω ∈ [0, π)) el periodograma superpuesto a la
% PSD del modelo que cumple con la ecuaci´on (6), en base a los par´ametros estimados

%genero ruido gaussiano blanco con N muestras

ruido = normrnd(0,1,1,N); 

%parametros para P = 2
[a_2, G_2] = ar_model(x, 2); %obtengo los coef G y a_i

%parametros para P = 13
[a_13, G_13] = ar_model(x, 13);

%parametros para P = 30
[a_30, G_30] = ar_model(x, 30);

%genero la señal de salida a partir de los parametros obtenidos para P={2, 12, 30}

%P=2
u_2 = G_2 * ruido;
x_2 = filter(G_2, [1 -a_2'], u_2); 

%P=13
u_13 = G_13 * ruido;
x_13 = filter(G_13, [1 -a_13'], u_13);

%P=30
u_30 = G_30 * ruido;
x_30 = filter(G_30, [1 -a_30'], u_30);


%calculo los periodogramas

%P=2
S_X2 = (1/N) * (abs(fft(x_2))).^2;

%P=12
S_X13 = (1/N) * (abs(fft(x_13))).^2;

S_X30 = (1/N) * (abs(fft(x_30))).^2;

%calculo la PSD teorica para cada señal

%P = 2

[H2, w] = freqz(G_2, [1 -a_2'], 'whole', N); 
SX2 = (abs(H2)).^2; %PSD teorica  

%P=13
H13 = freqz(G_13, [1 -a_13'], 'whole', N);
SX13 = (abs(H13)).^2; %PSD teorica 

%P=30
H30= freqz(G_30, [1 -a_30'], 'whole', N);
SX30 = (abs(H30)).^2; %PSD teorica 

%grafico los periodogramas
%P=2
figure();
plot(w/pi,10*log10(S_X2)); %psd para P=2
hold on 
plot(w/pi, 10*log10(SX2));
title("Densidad Espectral de Potencia");
legend("Periodograma con P=2", "PSD teorica con P=2", 'Location', 'north');
xlabel("$\omega/\pi$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');


%P=13
figure();
plot(w/pi,10*log10(S_X13)); %psd para P=13
hold on 
plot(w/pi, 10*log10(SX13));
title("Densidad Espectral de Potencia");
legend("Periodograma con P=13", "PSD teorica con P=13", 'Location', 'north');
xlabel("$\omega/\pi$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');

%P=30
figure();
plot(w/pi,10*log10(S_X30)); %psd para P=2
hold on 
plot(w/pi, 10*log10(SX30));
title("Densidad Espectral de Potencia");
legend("Periodograma con P=30", "PSD teorica con P=30", 'Location', 'north');
xlabel("$\omega/\pi$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');

%==========================================================================%
%Ejercicio 2
% Aplicando el m´etodo de Welch (Ap´endice 4.1), para una ventana de Hamming de largo
% M = 80 y un solapamiento del 50 %, estime la PSD de ambos registros de EEG. Grafique la
% PSD estimada con Welch junto a la PSD del modelo para P = 13.

% Parámetros del método de Welch
M = 80; % Ancho del segmento
overlap = M / 2; % Solapamiento del 50%


%estimamos la PSD para los datos del ECG con ojos abiertos mediante el
%metodo Welch

Pxx = metodo_Welch(x, N, M, overlap);

w1 = linspace(0,2, M);

% Graficar la PSD
figure();
plot(w1, 10*log10(Pxx));
hold on
plot(w/pi,10*log10(S_X13)); %psd para P=13
xlabel("$\omega/\pi$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
title('PSD usando el método de Welch');


% Utilice los par´ametros obtenidos del modelo para sintetizar las se˜nales EEG artificiales,
% suponiendo un orden P = 13. Luego estime los periodogramas mediante el m´etodo de Welch
% de forma analoga al ejercicio anterior. Compare en un gr´afico la PSD estimada para las se˜nales
% sintetizadas y la PSD del modelo te´orico seg´un la ecuaci´on (6).

%utilizamos la señal sintetizada con P=13 para estimar la PSD mediante el
%metodo Welch

Px13 = metodo_Welch(x_13', N, M, overlap);

figure();
plot(w1, 10*log10(Px13));
hold on
plot(w/pi, 10*log10(SX13)); %PSD teorica para P=13
xlabel("$\omega/\pi$",'Interpreter', 'latex');
ylabel("$S_X(\omega)$", 'Interpreter', 'latex');
title('PSD usando el método de Welch');




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
    Pxx_segments = zeros(M, L);

    for i = 1:L
        segment = x((i-1)*K + 1:(i-1)*K + M); % Seleccionar segmento
        segment = segment - mean(segment); % Remover la media del segmento
        window = hann(M); % Ventana de Hann
        segment = segment .* window; % Aplicar ventana
        Pxx_segments(:, i) = abs(fft(segment)).^2; % PSD del segmento
    end

    % Promediando las PSD de los segmentos
    Pxx = mean(Pxx_segments, 2);
    powV = mean(abs(window').^2);
    Pxx = Pxx/powV;

end
