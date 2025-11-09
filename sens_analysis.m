clear;
%% Question 5 - Sensitivity Analysis

% Defining time domain
dt = 0.01;
T = 40;     % Time duration of simulation
t = 0:dt:T;

% Constants
G = 30;        % Gait Force
f_init = 1;    % Gait Frequency (Hz)
N_init = 5000; % Number of People on the Bridge

m = 1e5;  % Bridge mass (kg)
c_init = 1.2e4;  % Bridge damping (kg/s)
k = 36e5; % Bridge stiffness (kg/s^2)

% Studied ranges for each variable
c_values = linspace(0,18e5,300);
N_values = linspace(1,5000,200);
f_values = linspace(0.7,1.2,200);

x_c = zeros(length(c_values),2*length(t)-1);
x_N = zeros(length(N_values),2*length(t)-1);
x_fmarch = zeros(length(f_values),2*length(t)-1);
x_fnorm = zeros(length(f_values),2*length(t)-1);

%% N Sensitivity Analysis
c = c_init;
f = f_init;

for Ni = 1:length(N_values)
    % Start the analysis with the updated parameters
    
    N = N_values(Ni);
    % Updated mass for the 2nd run, avg human weight of 60kg
    m = 70*N + 1e5;

    % Helpful Functions
    zeta = c/(2*sqrt(m*k));      % Damping Ratio (%)
    fn = (1/(2*pi)) * sqrt(k/m); % Natural Frequency (Hz)

    % Function F(s)
    Fs = 1;  % Initialising F(s)
    for i = 1:N
        % Fourier r(t) function for a unique individual

        % 1st Analysis
        % Produce a unique gait phase
        % gait_phase = rand() * 2 * pi;
        % rt = (4*G/pi)*(sin(2*pi*f*(t + gait_phase)) ...
        %        + (1/3)*sin(6*pi*f*(t + gait_phase)) ...
        %       + (1/5)*sin(10*pi*f*(t + gait_phase)) ...
        %       + (1/7)*sin(14*pi*f*(t + gait_phase)));

        % 2nd Analysis
        rt = (4*G/pi)*(sin(2*pi*f*(t)) ...
               + (1/3)*sin(6*pi*f*(t)) ...
              + (1/5)*sin(10*pi*f*(t)) ...
              + (1/7)*sin(14*pi*f*(t)));
        % Principle of superposition to sum up all the responses
        Fs = Fs + rt;
    end

    % Determinant calculation to determine the G(s) function
    determinant = c^2 - 4*m*k;

    if  abs(determinant) == 0
        alp = c/(2*m);
        Gs = t/m .* exp(-alp * t);                   % ∆ = 0
    elseif determinant > 0
        % Simplifying variables for clarity
        s1 = (-c + sqrt(determinant))/(2*m);
        s2 = (-c - sqrt(determinant))/(2*m);
        % G(s) equation
        Gs = (exp(s1*t) - exp(s2*t))./(m*(s1 - s2)); % ∆ > 0
    elseif determinant < 0
        % Simplifying variables for clarity
        alp = c/(2*m);
        bet = sqrt(-determinant)/(2*m);
        % G(s) equation
        Gs = (exp(-alp*t).*sin(bet*t))/(m*bet);      % ∆ < 0
    end

    % Calculated x function
    x_N(Ni,:) = conv(Fs, Gs, 'full') * dt;
end
figure;
maxN = max(x_N,[],2);
plot(N_values, maxN,'LineWidth',1.5, ...
    'DisplayName', sprintf('c = %.0f kg/s, f = %.2f Hz', c, f));
title('Sensitivity Analysis - Max Lateral Vibration against No. of People on Bridge');
legend('FontSize',14);
xlabel('Number of People');
ylabel('Maximum Lateral Displacement (m)');
grid on;

%% c Sensitivity Analysis
f = f_init;
N = N_init;

% Prepare the plot
for ci = 1:length(c_values)
    % Start the analysis with the updated parameters
    c = c_values(ci);

    % Helpful Functions
    zeta = c/(2*sqrt(m*k));      % Damping Ratio (%)
    fn = (1/(2*pi)) * sqrt(k/m); % Natural Frequency (Hz)
    
    % Function F(s)
    Fs = 0;  % Initialising F(s)
    for i = 1:N
        % Produce a unique gait phase
        gait_phase = rand() * 2 * pi;
        % Fourier r(t) function for a unique individual
        rt = (4*G/pi)*(sin(2*pi*f*(t+gait_phase)) ...
           + (1/3)*sin(6*pi*f*(t+gait_phase)) ...
           + (1/5)*sin(10*pi*f*(t+gait_phase)) ...
           + (1/7)*sin(14*pi*f*(t+gait_phase)));
        % Principle of superposition to sum up all the responses
        Fs = Fs + rt;
    end
    
    % Determinant calculation to determine the G(s) function
    determinant = c^2 - 4*m*k;
    
    if  abs(determinant) == 0
        alp = c/(2*m);
        Gs = t/m .* exp(-alp * t);                   % ∆ = 0
    elseif determinant > 0
        % Simplifying variables for clarity
        s1 = (-c + sqrt(determinant))/(2*m);
        s2 = (-c - sqrt(determinant))/(2*m);
        % G(s) equation
        Gs = (exp(s1*t) - exp(s2*t))./(m*(s1 - s2)); % ∆ > 0
    elseif determinant < 0
        % Simplifying variables for clarity
        alp = c/(2*m);
        bet = sqrt(-determinant)/(2*m);
        % G(s) equation
        Gs = (exp(-alp*t).*sin(bet*t))/(m*bet);      % ∆ < 0
    end
    
    % Calculated x function 
    x_c(ci,:) = conv(Fs, Gs, 'full') * dt;
end

% Plotting c Sensitivity Analysis
figure;
maxc = max(x_c,[],2);
plot(c_values, maxc,'LineWidth',1.5, ...
    'DisplayName', sprintf('N = %.0f, f = %.2f Hz', N, f));
title('Sensitivity Analysis - Max Lateral Vibration against damping coefficient');
legend('FontSize',14);
xlabel('Damping Coefficient (kg/s)');
ylabel('Maximum Lateral Displacement (m)');
grid on;

% Plotting c Sensitivity Analysis on a logarithmic scale
figure;
maxc = max(x_c,[],2);
semilogy(c_values, maxc,'LineWidth',1.5, ...
    'DisplayName', sprintf('N = %.0f, f = %.2f Hz', N, f));
title('Sensitivity Analysis - Max Lateral Vibration against damping coefficient');
legend('FontSize',14);
xlabel('Damping Coefficient (kg/s)');
ylabel('Maximum Lateral Displacement (m)');
grid on;

%% f Sensitivity Analysis
N = N_init;
c = c_init;

for fi = 1:length(f_values)
    % Start the analysis with the updated parameters
    f = f_values(fi);

    % Helpful Functions
    zeta = c/(2*sqrt(m*k));      % Damping Ratio (%)
    fn = (1/(2*pi)) * sqrt(k/m); % Natural Frequency (Hz)

    % Function F(s)
    Fs_march = 0;  % Initialising F(s)
    Fs_norm = 0;  % Initialising F(s)
    for i = 1:N
        % Produce a unique gait phase
        %gait_phase = rand() * 2 * pi;
        % Fourier r(t) function for a unique individual
        rt = (4*G/pi)*(sin(2*pi*f*(t)) ...
            + (1/3)*sin(6*pi*f*(t)) ...
            + (1/5)*sin(10*pi*f*(t)) ...
            + (1/7)*sin(14*pi*f*(t)));
        % Principle of superposition to sum up all the responses
        Fs_march = Fs_march + rt;

        gait_phase = rand() * 2 * pi;
        % Fourier r(t) function for a unique individual
        rt = (4*G/pi)*(sin(2*pi*f*(t+gait_phase)) ...
           + (1/3)*sin(6*pi*f*(t+gait_phase)) ...
           + (1/5)*sin(10*pi*f*(t+gait_phase)) ...
           + (1/7)*sin(14*pi*f*(t+gait_phase)));
        % Principle of superposition to sum up all the responses
        Fs_norm = Fs_norm + rt;
    end

    % Determinant calculation to determine the G(s) function
    determinant = c^2 - 4*m*k;

    if  abs(determinant) == 0
        alp = c/(2*m);
        Gs = t/m .* exp(-alp * t);                   % ∆ = 0
    elseif determinant > 0
        % Simplifying variables for clarity
        s1 = (-c + sqrt(determinant))/(2*m);
        s2 = (-c - sqrt(determinant))/(2*m);
        % G(s) equation
        Gs = (exp(s1*t) - exp(s2*t))./(m*(s1 - s2)); % ∆ > 0
    elseif determinant < 0
        % Simplifying variables for clarity
        alp = c/(2*m);
        bet = sqrt(-determinant)/(2*m);
        % G(s) equation
        Gs = (exp(-alp*t).*sin(bet*t))/(m*bet);      % ∆ < 0
    end

    % Calculated x function
    x_fmarch(fi,:) = conv(Fs_march, Gs, 'full') * dt;
    x_fnorm(fi,:) = conv(Fs_norm, Gs, 'full') * dt;
end
maxfmarch = max(x_fmarch,[],2);
maxfnorm = max(x_fnorm,[],2);

% Marching Figure
figure;
% Plot March
plot(f_values, maxfmarch,'-b','LineWidth',1.5, ...
    'DisplayName', sprintf('Marching Pace @ c = %.0f kg/s, N = %.0f', c, N));
% Plot natural frequency line
xline(fn, '--r', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Natural Frequency = %.4f Hz', fn));

legend('FontSize',14);
title('Sensitivity Analysis - Max Lateral Vibration against Gait Frequency');
xlabel('Gait Frequency (Hz)');
ylabel('Maximum Lateral Displacement (m)');
grid on;

% Normal Figure
figure;
% Plot Normal
plot(f_values, maxfnorm,'-b','LineWidth',1.5, ...
    'DisplayName', sprintf('Normal Pace @ c = %.0f kg/s, N = %.0f', c, N));
% Plot natural frequency line
xline(fn, '--r', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Natural Frequency = %.4f Hz', fn));

legend('FontSize',14);
title('Sensitivity Analysis - Max Lateral Vibration against Gait Frequency');
xlabel('Gait Frequency (Hz)');
ylabel('Maximum Lateral Displacement (m)');
grid on;
