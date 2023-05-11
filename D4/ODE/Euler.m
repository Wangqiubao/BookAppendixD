% Author-Zhang-Congqing
% Date-2023-4-20
clear;clc
y(1) = 1; %initial value
h = 0.1;  %step
t = 0:h:1;
% Euler
for i = 1:length(t)-1
    y(i+1)=y(i)+h*f(t(i),y(i));
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
