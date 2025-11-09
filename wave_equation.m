clear;
% Question 4 

% Defining time and space steps
dt = 0.1;
dx = 0.25;

% Defining time and space array axes
L = 1;
t = 0:dt:1;
x = 0:dx:L;

% Creating an empty mesh
u = zeros(length(x),length(t));

% Function u(i,j) in index form
u(1,:) = 0;         % u(0, t) = 0
u(end,:) = 0;       % u(L, t) = 0
u(:,1) = sin(pi*x); % u(x, 0) = f(x)

% First values determined by Neumann Boundary
for i = 2:length(x)-1
    u(i,2) = (21/25)*(sin(x(i)*pi)) + (2/25)*(u(i+1,1) + u(i-1,1));
end

% All other values infront of the known values
for j = 3:length(t)
    for i = 2:length(x)-1
      u(i,j) = (42/25)*u(i,j-1) + (4/25)*(u(i-1,j-1) ...
             + u(i+1,j-1)) - u(i,j-2);
    end
end

% Plotting the figure
figure;
surf(t,x,u);                    % 3D Surface Plot
set(gca, 'YDir', 'reverse');    % Flip x axis to go from 0 to 1
title("Function u against displacement and time");
xlabel('Time, t (s)');
ylabel('Length, x (m)');
zlabel('Amplitude');