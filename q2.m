clear;
% Question 2

% Creating time domain
T = 3;
dt = 0.01;
t = 0:dt:T;

% Other Constants
f = 1;             % Gait Frequency
G = 30;            % Gait Force Magnitude

% Square wave function
s = G*(2*(mod(t,1/f) < (1/f)/2) - 1);

% Fourier series function
r = (4*G/pi)*(sin(2*pi*t*f) ...
    + (1/3) *sin(6*pi*t*f) ...
    + (1/5)*sin(10*pi*t*f) ...
    + (1/7)*sin(14*pi*t*f));

% Plot the two graphs on the same figure
figure;
plot(t, r, 'LineWidth', 1.5, ...
    'DisplayName', 'First 4 non-zero terms Fourier Expansion');
hold on;
plot(t, s, 'LineWidth', 1.5, 'DisplayName', 'Square Wave');
hold off;
xlabel('Time (s)');
ylabel('Lateral Walking Force (N)');
ylim([-40 55]);
title('Fourier Series Expansion of Periodic Force Function, r(t)');
legend('FontSize',12);
grid on;