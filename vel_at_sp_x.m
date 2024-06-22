clc;
clear all;

beta = 0.01;
lamda = pi / 2;
b = 4 / pi;
mobility = -10^(-10);
Ds = 10^(-9);
L = 10^(-4);


X_positions = [0.01, 0.5, 1];

time = linspace(0, 1, 10);  

v_profile_01 = zeros(1, length(time));
v_profile_05 = zeros(1, length(time));
v_profile_1 = zeros(1, length(time));

for n = 1:length(time)
    t = time(n);
    T = t * Ds / L^2;
    
    v_profile_01(n) = (mobility * lamda * b * (1 - beta) * cos(lamda * X_positions(1)) * exp(-lamda^2 * T)) / ...
                      (Ds * (beta + (1 - beta) * sin(lamda * X_positions(1)) * exp(-lamda^2 * T)));
                  
    v_profile_05(n) = (mobility * lamda * b * (1 - beta) * cos(lamda * X_positions(2)) * exp(-lamda^2 * T)) / ...
                      (Ds * (beta + (1 - beta) * sin(lamda * X_positions(2)) * exp(-lamda^2 * T)));
                  
    v_profile_1(n) = (mobility * lamda * b * (1 - beta) * cos(lamda * X_positions(3)) * exp(-lamda^2 * T)) / ...
                     (Ds * (beta + (1 - beta) * sin(lamda * X_positions(3)) * exp(-lamda^2 * T)));
end

figure;
plot(time, v_profile_01, '-o', 'DisplayName', ['x = ', num2str(X_positions(1))]);
hold on;
plot(time, v_profile_05, '-x', 'DisplayName', ['x = ', num2str(X_positions(2))]);
plot(time, v_profile_1, '-s', 'DisplayName', ['x = ', num2str(X_positions(3))]);
hold off;

legend;
xlabel('Time');
ylabel('Velocity v');
title('Velocity Profile over Time at Specific Positions');
grid on;