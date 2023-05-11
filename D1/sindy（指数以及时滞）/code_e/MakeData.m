function [X,t,Delay] = MakeData()
%%
% f 和 g  参数
% 根据自己的方程进行设置
% 若无时滞项 则 tau = 0

tau = 0;
m = 1;

%%
% h 为步长
% N 为计算次数
% D 为噪声强度矩阵 
% D(1) 为GSN1的噪声强度，D(2) 为GSN2的噪声强度
% f0、g0分别为 x1 和 x2的初值
f0 = 1;
g0 = -1;
h = 0.001;
N = 10001;
D = [1,1];
Delay = tau/h;
%%
% 设置方程
% f -- x1'
% g -- x2'
f = @(x1tau,GSN1,GSN2,x1,x2)(x2 + GSN2);
g = @(x1tau,GSN1,GSN2,x1,x2)((5 - exp(x1) + x1^2)/m + GSN1 );
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