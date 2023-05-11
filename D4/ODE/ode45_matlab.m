% Author-Zhang-Congqing
% Date-2023-4-20
clear;clc
y0 = 1; %initial value
tspan = [0,1];
% Runge-Kutta
[t,y] = ode45(@f,tspan,y0);
%% plot
plot(t,exp(-5*t))
hold on
plot(t,y,'k*')
legend('Exact solution','Euler')
%% ODE
function dy = f(t,y)
    dy = -5*y;
end
