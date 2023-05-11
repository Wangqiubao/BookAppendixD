
clear all; close all; clc
figpath = '../figures/';
addpath('./utils');

%% generate Data
polyorder = 2;  % search space up to fifth order polynomials
usesine = 0;    % no trig functions
usee = 1;
n = 2;          % 2D system
[x,tspan,delay] = MakeData();
%%
x0 = [1, -1];        % 初值，保持和MakeData中一致
%% compute Derivative 
eps = 0.1;      % noise strength
dx = [gradient(x(:,1))./gradient(tspan'),gradient(x(:,2))./gradient(tspan')];
dx = dx + eps*randn(size(dx));   % add noise
figure(1)
plot(tspan,x(:,1));
hold on
plot(tspan,x(:,2))
hold on
%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine,delay,usee);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambda = 0.5;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n)
poolDataLIST({'x','y'},Xi,n,polyorder,usesine,delay,usee);
syms x1 x2
yin = [x1,x2];
%% integrate true and identified systems
[f,g] = sparseGalerkin(yin,Xi,polyorder,usesine,delay,usee);
h = 0.001;
N = 10001;
D = [1,1];
tau = 0;
[tx1,tx2] = Runge_Kutta(h,N,f,g,tau,D,x0(1),x0(2)); 
figure(2)
plot(x(:,1),x(:,2),'r-');
hold on
plot(tx1,tx2,'k--');

figure(1)
plot(tspan,tx1,'k--')
hold on
plot(tspan,tx2,'k--')
legend('x1','x2','拟合曲线')