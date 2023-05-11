function [X,t,Delay] = MakeData()
%%
% f 和 g  参数
% 根据自己的方程进行设置
% 若无时滞项 则 tau = 0

p = 0.1;
T = 100;
tau = 3;
eta = 0.5;
epsilon = 0.27;
alpha = sqrt((1-eta)/(2*p*T));
beta = 2*p/(1-eta);

%%
% h 为步长
% N 为计算次数
% D 为噪声强度矩阵 
% D(1) 为GSN1的噪声强度，D(2) 为GSN2的噪声强度
% f0、g0分别为 x1 和 x2的初值
f0 = 0.1;
g0 = 0;
h = 0.001;
N = 20001;
D = [1,1];
Delay = tau/h;
%%
% 设置方程
% f -- x1'
% g -- x2'
f = @(x1tau,GSN1,GSN2,x1,x2)(x2 + epsilon*x1*x2);
g = @(x1tau,GSN1,GSN2,x1,x2)(-x1 - alpha*(1+beta)*x2 + eta*x1tau - epsilon*alpha*beta*x1*x2);
[X1,X2] = Runge_Kutta(h,N,f,g,tau,D,f0,g0); 
X = [X1;X2]';
%%
t = 0:h:(N-1)*h;
% % 绘图
% figure(1)
% plot(t,X1);
% xlabel t;ylabel x;
% figure(2)
% plot(t,X2);
% xlabel t;ylabel y;
% 
% figure(3)
% plot(X1,X2);
% xlabel x;ylabel y;