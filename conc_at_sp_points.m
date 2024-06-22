clc;
clear all;


beta = 0.01;
co = 1;
c1 = beta;
b = 4 / pi;
Ds = 10^(-9); 
L = 10^(-4);
t = 10;

X = linspace(0, 1, 10);


T1 = .4;
T2 = .5 ;
T3 = 1;


c1_profile = zeros(1, length(X));
c2_profile = zeros(1, length(X));
c3_profile = zeros(1, length(X));


for m = 1:length(X)
    c1_profile(m) = beta + (1-beta)*b*sin(pi*X(m)/2)*exp(-pi^2*T1*0.25);
    c2_profile(m) = beta + (1-beta)*b*sin(pi*X(m)/2)*exp(-pi^2*T2*0.25);
    c3_profile(m) = beta + (1-beta)*b*sin(pi*X(m)/2)*exp(-pi^2*T3*0.25);
end


figure;
plot(X, c1_profile, '-o',LineWidth=2);
hold on;
plot(X, c2_profile, '-x',LineWidth=2);
hold on
plot(X, c3_profile, '-s',LineWidth=2 );
legend('at 0.4td','at 0.5td','at td',Location='best')
xlabel('X');
ylabel('Concentration c');
title('Concentration Profile at Specific Times');