clear all, close all, clc
addpath('utils\')
% load conicity.mat;
% load conicity_ploy.mat;

polyorder = 1;
n = 12;
usesine = 0;
dt = 0.01;
tspan=dt:dt:100;
options = odeset('RelTol',1e-8,'AbsTol',1e-8*ones(1));

y = [0.001,0,0.001,0,0.001,0,0.001,0,0.001,0,0.001,0];
[t,x]=ode45(@(t,x) fun(t,x,98),tspan,y,options);

%% compute Derivative
eps = 1;
% for i = 1:length(y)
%     dx(:,i) = gradient(x(:,i))./gradient(tspan');
% end

for i=1:length(x)
    dxtmp(i,:) = fun(0,x(i,:),98);
end
dx = dxtmp(1:1000,:);
%% pool Data  (i.e., build library of nonlinear time series)
Theta = poolData(x(1:1000,:),n,polyorder,usesine);
m = size(Theta,2);

%% compute Sparse regression: sequential least squares
lambda0 = 0.1;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda0,n);

%% FIGURE 1:  LORENZ for T\in[0,20]
tspan = [0 20];
[tA,xA]=ode45(@(t,x)fun(t,x,98),tspan,y,options);   % true model
[tB,xB]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polyorder,usesine),tspan,y,options);  % approximate
figure(1)
plot(tA(end-1000:end,1),xA(end-1000:end,1),'k-',tB(end-1000:end,1),xB(end-1000:end,1),'r.')
figure(2)
plot(xA(end-1000:end,1),xA(end-1000:end,5),'k-');
hold on
plot(xB(end-1000:end,1),xB(end-1000:end,5),'r.');

function dx = fun(t,y,V)
    dy1 = y(2);
    dy2 = 6944.444444*y(11) + (50000*y(9))/9 - (50000*y(1))/9 - (140000*y(2))/(9*V) + (140000*y(5))/9 - 81.76676339*lambda(y(1))*y(1) - Ft(y(1))/1800;
    
    dy3 = y(4);
    dy4 = 6944.444444*y(11) + (50000*y(9))/9 - (50000*y(3))/9 - (140000*y(4))/(9*V) + (140000*y(7))/9 - 81.76676339*lambda(y(3))*y(3) - Ft(y(3))/1800;
    
    dy5 = y(6);
    dy6 = -178354.2857*y(5) + 178354.2857*y(11) - 69441.86049*lambda(y(1))*y(1) - 22290.49000*y(6)/V + 117.1685071*lambda(y(1))*y(5);
    
    dy7 = y(8);
    dy8 = -178354.2857*y(7) + 178354.2857*y(11) - 69441.86049*lambda(y(3))*y(3) - 22290.49000*y(8)/V + 117.1685071*lambda(y(3))*y(7);
    
    dy9 = y(10);
    dy10 = -(300*y(10))/23 - 8767.826087*y(9) + (100000*y(1))/23 + (100000*y(3))/23;
    
    dy11 = y(12);
    dy12 = -187.5721*y(12) - 108113.7750*y(11) + 48018.46154*y(5) + 48018.46154*y(7) + 4807.692308*y(1) - 4807.692308*y(3);
    dx = [dy1;dy2;dy3;dy4;dy5;dy6;dy7;dy8;dy9;dy10;dy11;dy12;];
end

function F_t = Ft(x)
    K_r = 1.617 * 10^7;
    if  x > 0.00923
        F_t = K_r*(x - 0.00923);
    elseif (x >= -0.00923) & (x <= 0.00923)
        F_t = 0;
    else
        F_t = K_r*(x + 0.00923);
    end
    
end


function conicity = lambda(x)
%     x = x * 1000;
%     con_poly = [-6.50477366354501e-05,0.00177006663721497,-0.0194568730625890,0.112742689482152,-0.368659901106896,0.639057225082614,-0.276497119936299];
%     conicity = con_poly * [x^6;x^5;x^4;x^3;x^2;x;1];
%     if conicity < 0.1
%         conicity = 0.104;
%     end
    conicity = 0.2;
end 