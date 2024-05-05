clear all
clc

pkg load optim
pkg load miscellaneous
pkg load statistics
pkg load signal

% Generar la señal
ojosAbiertos = dlmread('eeg_ojos_abiertos_t7.csv','\t',0,0);
x = ojosAbiertos';
N = length(x); % Número de muestras
fs = 200; % Frecuencia de muestreo
t = 0:1/fs:N/fs; % Tiempo

% Parámetros del método de Welch
M = 80; % Ancho del segmento
overlap = M / 2; % Solapamiento del 50%

% Calcular PSD usando el método de Welch
K = M - overlap; % Distancia entre segmentos
    
% Calculando la PSD usando el método de Welch
L = floor((N - M) / K) + 1; % Número total de segmentos
Pxx_segments = zeros(M, L);
    
for i = 1:L
    segment = x((i-1)*K + 1:(i-1)*K + M); % Seleccionar segmento
    segment = segment - mean(segment); % Remover la media del segmento
    window = hann(M); % Ventana de Hann
    segment = segment .* window'; % Aplicar ventana
    Pxx_segments(:, i) = abs(fft(segment)).^2; % PSD del segmento
end
    
% Promediando las PSD de los segmentos
Pxx = mean(Pxx_segments, 2);
powV = mean(abs(window').^2);
Pxx = Pxx/powV;
    
% Calculando las frecuencias correspondientes
f = (0:M-1) * fs / M;

% Graficar la PSD
plot(f, 10*log10(Pxx));
xlabel('Frecuencia (Hz)');
ylabel('Densidad Espectral de Potencia (dB/Hz)');
title('PSD usando el método de Welch');
