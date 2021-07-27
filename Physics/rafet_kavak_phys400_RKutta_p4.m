% Rafet KAVAK - 2166783
% PHYS400 Project - 22.06.2021
% v2 - 12.07.2021

clc; % Clears the Command Window
close all; % Closes the all figures
clear; % Clears the Workspace

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%
tStart = tic;

% f and Initial Conditions
f1 = @(t,x)(x(3));
f2 = @(t,x)(x(4));
f3 = @(t,x)(-x(1)*(x(6)+x(2)^2)-2*x(7)*x(2));
f4 = @(t,x)(-x(2)*(x(5)+x(1)^2)-2*x(7)*x(1));
f5 = @(t,x)(2*x(8));
f6 = @(t,x)(2*x(9));
f7 = @(t,x)(x(11)+x(10));
f8 = @(t,x)(x(12)-x(5)*(x(6)+x(2)^2)-2*x(7)*(x(7)+x(1)*x(2)));
f9 = @(t,x)(x(13)-x(6)*(x(5)+x(1)^2)-2*x(7)*(x(7)+x(1)*x(2)));
f10 = @(t,x)(x(14)-x(7)*(x(5)+x(1)^2)-2*x(5)*(x(7)+x(1)*x(2)));
f11 = @(t,x)(x(14)-x(7)*(x(6)+x(2)^2)-2*x(6)*(x(7)+x(1)*x(2)));
f12 = @(t,x)(-2*x(8)*(x(6)+x(2)^2)-4*x(11)*(x(7)+x(1)*x(2)));
f13 = @(t,x)(-2*x(9)*(x(5)+x(1)^2)-4*x(10)*(x(7)+x(1)*x(2)));
f14 = @(t,x)(-x(10)*(x(6)+x(2)^2)-2*x(9)*(x(7)+x(1)*x(2))-x(11)*(x(5)+x(1)^2)-2*x(8)*(x(7)+x(1)*x(2)));

f  = @(t,x)([f1(t,x);f2(t,x);f3(t,x);f4(t,x);f5(t,x);f6(t,x);f7(t,x);f8(t,x);f9(t,x);f10(t,x);f11(t,x);f12(t,x);f13(t,x);f14(t,x)]);

factor = 4;

x0 = [0.625*factor; %x1
      0.325*factor;
        0;
        0;
       0.5;
       0.5;
        0;
        0;
        0;
        0;
        0;
       0.5;
       0.5;
        0];
    

% Sim Parameters
tend  = 6;
h     = 0.001;      % initial step size

% Simulation
clear x time;
x(:,1) = x0;                % state values during the simulation
time(1) = 0;                % simulation time

for kk = 1:tend/h
    x(:,kk+1) = RKmP4(time(kk),x(:,kk),h,f); 
    time(kk+1) = kk*h;
end

figure; 
plot(time,x(1,:),'LineWidth',2), grid; 
xlabel('{\boldmath$t/\tau_{L}$}'), ylabel('{\boldmath$<\hat{x}(t)>$}'); 
factorstr = sprintf('$f = %.0f$', factor);
legend(factorstr); 
title('$dx/dt = f(x)$');

%for f = 1: ylim([-0.8 0.8]);
%for f = 2: ylim([-1.5 1.5]);
%for f = 4: ylim([-2 3]);
%for f = 8: ylim([-3 6]);
ylim([-3 6]);

tEnd = toc(tStart)

function xNext = RKmP4(t,x,h,f) 
c2 = 1/2; 
c3 = 1/2; 
c4 = 1; 
b1 = 1/6; 
b2 = 2/6; 
b3 = 2/6; 
b4 = 1/6; 
a21 = 1/2; 
a31 = 0; 
a32 = 1/2; 
a41 = 0; 
a42 = 0; 
a43 = 1; 

k1 = f(t, x); 
k2 = f(t+c2*h, x+h*a21*k1); 
k3 = f(t+c3*h, x+h*a31*k1+h*a32*k2); 
k4 = f(t+c4*h, x+h*a41*k1+h*a42*k2+h*a43*k3); 
xNext = x+h*(b1*k1+b2*k2+b3*k3+b4*k4); 
end


