
Ds=10^(-9);
Dp=10^(-12);
mob=-10^(-10);
l=[pi/2 3*pi/2 5*pi/2];
L=10^(-4);
beta=0.01;
b=2./l;

t=linspace(0,10,10001);
x=[0.2  0.5 1];
t_scale1= exp(-1.*(l(1).^2).*t);
t_scale2= exp(-1.*(l(2).^2).*t);
 t_scale3=exp(-1.*(l(3).^2).*t);
 dcdx=(l(1).*(1-beta).*b(1)*cos(l(1).*x(1)).*t_scale1 );
vp_1=mob.*(l(1).*(1-beta).*b(1)*cos(l(1).*x(1)).*t_scale1 +l(2).*(1-beta).*b(2)*cos(l(2).*x(1)).*t_scale2 +l(3).*(1-beta).*b(3)*cos(l(3).*x(1)).*t_scale3)./(Ds*(beta +(1-beta).*(b(1)*sin(l(1).*x(1)).*t_scale1 +b(2).*sin(l(2).*x(1)).*t_scale2+b(3)*sin(l(3).*x(1)).*t_scale3)));
vp_1FTA=mob.*(l(1).*(1-beta).*b(1)*cos(l(1).*x(1)).*t_scale1)./(Ds*(beta +(1-beta).*(b(1)*sin(l(1).*x(1)).*t_scale1)));


figure
 plot(t,vp_1,'k',linewidth=2.5);
hold on
plot(t,vp_1FTA,'--r',linewidth=2.5);
xlim([0,0.2]);

legend ('at x=0.2','at x=0.2FTA',Location='best');%CHANGE THE TIME SCALE TITLE
xlabel('characteristic time');
ylabel('velocity')
grid minor
title('Velocity profile curve');