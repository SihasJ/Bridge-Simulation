clear;
% Question 1 

% Creating time domain
T = 4;
dt = 0.01;
t = 0:dt:T;

% Defining transform functions
F = sin(2 * pi * t);
G = t .* exp(-6 * t);
H = sin(10 * pi * t);

% Calculated x function with dt scaling correction
x1 = 1e-5 * conv(F, G, 'full') * dt;
x2 = 1e-5 * conv(H, G, 'full') * dt;

% Matching convolution result with new time domain
t_conv = linspace(0, 2*t(end), length(x1));

% Plot 1
figure;
plot(t_conv, x1,'b-','LineWidth',1,'DisplayName','$10^{-5} \left( \sin(2 \pi t) t e^{-6t} \right)$');
title('Lateral Vibrations of a Bridge');
xlabel('Time (s)');
ylabel('Displacement from centre (m)');
xlim([0 T]);
legend('Interpreter','latex','FontSize',14);
grid on;

% Plot 2
figure;
plot(t_conv, x2,'r-','LineWidth',1,'DisplayName','$10^{-5} \left( \sin(10 \pi t) t e^{-6t} \right)$');
title('Lateral Vibrations of a Bridge');
xlabel('Time (s)');
ylabel('Displacement from centre (m)');
xlim([0 T]);
legend('Interpreter','latex','FontSize',14);
grid on;

% Plot 3
figure;
plot(t_conv, x1,'b-','LineWidth',1,'DisplayName','$10^{-5} \left( \sin(2 \pi t) t e^{-6t} \right)$');
hold on;
plot(t_conv, x2,'r-','LineWidth',1,'DisplayName','$10^{-5} \left( \sin(10 \pi t) t e^{-6t} \right)$');
hold off
title('Lateral Vibrations of a Bridge');
xlabel('Time (s)');
ylabel('Displacement from centre (m)');
xlim([0 T]);
legend('Interpreter','latex','FontSize',14);
grid on;