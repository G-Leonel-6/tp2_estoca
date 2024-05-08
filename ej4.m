rp = 0.01;
rs = 0.001;
M = 10000;
fs = 200;

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

figure
hold on
grid on
title('Filtro de ondas Delta')
xlim([0, 45])
line([0 3], [1.01 1.01], color = 'red', linestyle = '--')
line([0 3], [0.99 0.99], color = 'red', linestyle = '--')
line([5 10], [0.001 0.001], color = 'red', linestyle = '--')
line([5 10], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{D}(f)|')
plot(w, abs(H), color = 'blue')
hold off

figure
hold on
grid on
title('Filtro de ondas Theta')
xlim([0, 45])
line([5 8], [1.01 1.01], color = 'red', linestyle = '--')
line([5 8], [0.99 0.99], color = 'red', linestyle = '--')
line([0 3], [0.001 0.001], color = 'red', linestyle = '--')
line([0 3], [-0.001 -0.001], color = 'red', linestyle = '--')
line([10 12], [0.001 0.001], color = 'red', linestyle = '--')
line([10 12], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{T}(f)|')
plot(w, abs(H1), color = 'blue')
hold off

figure
hold on
grid on
title('Filtro de ondas Alpha')
xlim([0, 45])
line([10 13], [1.01 1.01], color = 'red', linestyle = '--')
line([10 13], [0.99 0.99], color = 'red', linestyle = '--')
line([0 8], [0.001 0.001], color = 'red', linestyle = '--')
line([0 8], [-0.001 -0.001], color = 'red', linestyle = '--')
line([15 20], [0.001 0.001], color = 'red', linestyle = '--')
line([15 20], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{A}(f)|')
plot(w, abs(H2), color = 'blue')
hold off

figure
hold on
grid on
title('Filtro de ondas Beta')
xlim([0, 45])
line([15 29], [1.01 1.01], color = 'red', linestyle = '--')
line([15 29], [0.99 0.99], color = 'red', linestyle = '--')
line([0 13], [0.001 0.001], color = 'red', linestyle = '--')
line([0 13], [-0.001 -0.001], color = 'red', linestyle = '--')
line([31 35], [0.001 0.001], color = 'red', linestyle = '--')
line([31 35], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{B}(f)|')
plot(w, abs(H3), color = 'blue')
hold off

figure
hold on
grid on
title('Filtro de ondas Gamma')
xlim([0, 45])
plot(w, abs(H4), color = 'blue')

line([31 45], [1.01 1.01], color = 'red', linestyle = '--')
line([31 45], [0.99 0.99], color = 'red', linestyle = '--')
line([0 29], [0.001 0.001], color = 'red', linestyle = '--')
line([0 29], [-0.001 -0.001], color = 'red', linestyle = '--')
xlabel('f[Hz]')
ylabel('|H_{G}(f)|')

hold off