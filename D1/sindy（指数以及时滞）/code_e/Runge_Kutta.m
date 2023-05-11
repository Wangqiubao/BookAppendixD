function [X1,X2] = Runge_Kutta(h,N,f,g,tau,D,f0,g0)
%%
%                                    %
% 计算带有时滞/随机/碰撞的微分方程组   %
%                                    %
%%
% h——步长
% N——计算次数
% f,g——待计算函数
% tau——时滞
% D——噪声强度
% f0、g0——初值
%%
%其它参数设置
% r = 0.95;%碰撞能量损失系数，默认0.95

%生成高斯白噪声
GSN = [];
for i = 1:length(D)
    GsnTemp = wgn(1,N,D(i)) * sign(D(i));
    GSN = [GSN;GsnTemp]; 
end
%%
%设置初值
X1 = ones(1,N)*f0;
X2 = ones(1,N)*g0;
%%
% Runge-Kutta 4
for i = 1:1:N-1
    %时滞
    if i*h <= tau
        x1tau = 0;
    else
        x1tau = X1(ceil(i-(tau-mod(tau,h))/h));
    end
    x1 = X1(i);x2 = X2(i);
    Kf1 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1,x2);
    Kg1 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1,x2);
    Kf2 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1+Kf1*h/2,x2+Kg1*h/2);
    Kg2 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1+Kf1*h/2,x2+Kg1*h/2);
    Kf3 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1+Kf2*h/2,x2+Kg2*h/2);
    Kg3 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1+Kf2*h/2,x2+Kg2*h/2);
    Kf4 = feval(f,x1tau,GSN(1,i),GSN(2,i),x1+Kf3*h,x2+Kg3*h);
    Kg4 = feval(g,x1tau,GSN(1,i),GSN(2,i),x1+Kf3*h,x2+Kg3*h);

    X1(i+1) = X1(i) + (Kf1 + 2*Kf2 + 2*Kf3 + Kf4)*h / 6;
    X2(i+1) = X2(i) + (Kg1 + 2*Kg2 + 2*Kg3 + Kg4)*h / 6;
end