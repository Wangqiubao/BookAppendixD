
clear all; close all; clc
figpath = '../figures/';
addpath('./utils');

%% generate Data
polyorder = 2;  % search space up to fifth order polynomials
usesine = 0;    % no trig functions
n = 2;          % 2D system
[x,tspan,delay] = MakeData();
%%
x0 = [0.1, 0];        % 初值，保持和MakeData中一致
options = odeset('RelTol',1e-10,'AbsTol',1e-10*ones(1,n));
%% compute Derivative 
eps = 0.001;      % noise strength
dx = [gradient(x(:,1))./gradient(tspan'),gradient(x(:,2))./gradient(tspan')];
dx = dx + eps*randn(size(dx));   % add noise
figure()
plot(dx(:,1));
hold on
plot(dx(:,2))
usetau = 1;
%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x,n,polyorder,usesine,delay);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambda = 0.015;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n)
poolDataLIST({'x','y'},Xi,n,polyorder,usesine,delay);
syms x1 x2
yin = [x1,x2];
%% integrate true and identified systems
[f,g] = sparseGalerkin(yin,Xi,polyorder,usesine,usetau);
f0 = 0.1;
g0 = 0;
h = 0.001;
N = 20001;
D = [1,1];
tau = 3;
[tx1,tx2] = Runge_Kutta(h,N,f,g,tau,D,f0,g0); 
figure()
plot(x(:,1),x(:,2),'r-');
hold on
plot(tx1,tx2,'k--');