% DDE No boundary restrictions
% Author-Zhang-Congqing
% Date-2023-4-20
clear;clc;close all;
timebgn = 0;
timeend = 100;
tspan = [timebgn,timeend];   % time interval
tau = [0.1 1];               % time delay
y0 = [0.1;1];                % initial value
sol=dde23(@myddefun, tau, y0, tspan);

% plot
subplot(221)
plot(sol.x,sol.y(1,:))
subplot(223)
plot(sol.x,sol.y(2,:))
subplot(122)
plot(sol.y(1,:),sol.y(2,:))
%% SDDE setting
function dy = myddefun(t,y,Z) % equation being solved
    %The parameters in the equation are written here
    delta = 0.2;w0 = 1;alpha = 1;beta = 2;epsilon = 0.1;
    ylag1 = Z(:,1);
    ylag2 = Z(:,2);
    dy = zeros(2,1);           % a column vector
    dy(1) = y(2);
    dy(2) = (beta-w0^2)*y(1)-delta*y(2)- beta/2*(ylag1(1)+ylag2(1)) -epsilon*alpha*y(1)^3;
end

