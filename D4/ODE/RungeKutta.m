% Author-Zhang-Congqing
% Date-2023-4-20
clear;clc
y(1) = 1; %initial value
h = 0.1;  %step
t = 0:h:1;
% Runge-Kutta
for i = 1:length(t)-1
    k1 = f(t(i),y(i));
    k2 = f(t(i)+h/2,y(i)+h*k1/2);
    k3 = f(t(i)+h/2,y(i)+h*k2/2);
    k4 = f(t(i)+h,y(i)+h*k3);
    y(i+1)=y(i)+h*(k1+2*k2+2*k3+k4)/6;
end
%% plot
plot((0:h/10:1),exp(-5*(0:h/10:1)))
hold on
plot(t,y,'k*-')
legend('Exact solution','Euler')
%% ODE
function dy = f(t,y)
    dy = -5*y;
end
