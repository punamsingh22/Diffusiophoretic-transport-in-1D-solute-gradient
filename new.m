clc;
clear all;

M = 51; %for number of space nodes
N = 51; %for number of time nodes
mob = -10^(-10); %non-zero mobility
% mob = 0        %zero mobility
Ds = 10^(-9);
Dp = 10^(-12);
beta = 0.01;
lamda = pi/2;

co = 1;
c1 = beta;
b = 4/pi;

X = linspace(0, 1, M);
T = linspace(0, 1, N);

c = zeros(length(T), length(X));

for t = 1:length(T)
    for m = 1:length(X)
        c(t, m) = beta + (1 - beta) * b * sin(pi * X(m) / 2) * exp(-(pi)^2 * T(t) * 0.25);
    end
end

v = zeros(length(T), length(X));


for t = 1:length(T)
    for m = 1:length(X)
        v(t,m) = (mob*lamda*b*(1-beta)*cos(lamda*X(m))*(exp(-lamda*lamda*T(t))))/(Ds*(beta+(1-beta)*sin(lamda*X(m))*exp(-lamda*lamda*T(t))));   
    end
end

dx = X(2);
dt = T(2);

n0 = zeros(length(T), length(X));

fun = @(n) root (n,M,N,v,T,X);
n = fsolve(fun, n0);
figure(1);
hold on;
plot(X, n(floor(N * 0.5), :),LineWidth=2);  % Plot at T = 0.5
plot(X, n(floor(N * 0.4), :),LineWidth=2);  % Plot at T = 0.4
plot(X, n(floor(N * 0.3), :),LineWidth=2);  % Plot at T = 0.3
plot(X, n(floor(N * 0.2), :),LineWidth=2);  % Plot at T = 0.2
xlabel('X');
ylabel('Particle concentration n');
title('Particle concentration profiles at different times');
legend('T = 0.5', 'T = 0.4', 'T = 0.3', 'T = 0.2',Location='best');
hold off;

% figure;
% hold on;
% plot(X, n(floor(N * 5/5), :)); % Plot at T = 5
% plot(X, n(floor(N * 4/5), :)); % Plot at T = 4
% plot(X, n(floor(N * 3/5), :)); % Plot at T = 3
% plot(X, n(floor(N * 2/5), :)); % Plot at T = 2
% xlabel('X');
% ylabel('Particle concentration n');
% title('Particle concentration profiles at different times');
% legend('T = 5', 'T = 4','T = 3','T = 2');
% hold off;

function eq = root(n,M,N,v,T,X)
    no = 0;
    n1 = 1;
  
    Ds = 10^(-9);
    Dp = 10^(-12);
    
    dx = X(2);
    dt = T(2);
    
    eq = zeros(length(T), length(X));  
    
    
    eq(1,:) = n(1,:) - n1;
    eq(:,1) = n(:,1) - no;
    for t = 1:N-1
        for m = 2:M-1
            eq(t+1,m) = (n(t+1,m)-n(t,m))/dt + v(t+1,m)*(n(t+1,m+1)-n(t+1,m))/dx + n(t+1,m)*(v(t+1,m+1)-v(t+1,m))/dx - (Dp/Ds)*(n(t+1,m-1)-2*n(t+1,m)+n(t+1,m+1))/dx^2;
        end
        eq(t+1,M) = n(t+1,M) - n(t+1,M-1);
    end
end



