clc;
clear all;

M = 10;
L = 10^(-4);
mobility = -10^(-10);
Ds = 10^(-9);
Dp = 10^(-12);
beta = 0.01;
lamda = pi/2;
b = 4/pi;

Time = 0.050;
dx = L/(M-1);
dt = Time/10;
I = Time/dt+1;


n0 = zeros(I, M);

fun = @root;
n = fsolve(fun, n0);


function eq = root(n)
    M = 10;
    no = 0;
    n1 = 1;
    L = 10^(-4);
    mobility = -10^(-10);
    Ds = 10^(-9);
    Dp = 10^(-12);
    beta = 0.01;
    lamda = pi/2;
    b= 4/pi;
    Time = 0.050;
    dx = L/(M-1);
    dt = Time/10;
    I = Time/dt+1;

    eq = zeros(I, M);  
    
    eq(:,1) = n(:,1) - no;
    eq(1,:) = n(1,:) - n1;

    for t = 1:I-1
        for m = 2:M-1
            eq(t+1,m) = (n(t+1,m)-n(t,m))/dt...
                + (mobility*lamda*b*(1-beta)*cos(lamda*dx*(m-1))*(exp(-lamda*lamda*dt*(t-1))*(n(t+1,m+1)-n(t+1,m))))/(Ds*dx*(beta+(1-beta)*sin(lamda*dx*(m-1))*exp(-lamda*lamda*dt*(t-1)))) ...
               - n(t+1,m)*(mobility*lamda*lamda*b*(1-beta)*exp(-1*lamda*lamda*dt*(t-1))*(beta*sin(lamda*dx*(m-1))+(1-beta)*exp(-lamda*lamda*dt*(t-1))))/(Ds*(beta+(1-beta)*sin(lamda*dx*(m-1))*exp(-lamda*lamda*dt*(m-1)))^2)...
                - (Dp*(n(t+1,m-1)-2*n(t+1,m)+n(t+1,m+1)))/(Ds*(dx)^2);
        end
        eq(t+1,M) = n(t+1,M) - n(t+1,M-1);
    end
end


figure(1);
hold on;
legends = cell(1, I);
for i = 1:I
    plot(n(i,:));
    legends{i} = ['Time Step ', num2str(i)];
end
legend(legends);
hold off;
