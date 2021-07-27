% Rafet KAVAK - 2166783
% PHYS400 Project - 22.06.2021
% v2 - 12.07.2021

clc; % Clears the Command Window
close all; % Closes the all figures
clear; % Clears the Workspace

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

%{
Variables
------------------------------
    x      ------> x1
    y      ------> x2
   p_x     ------> x3
   p_y     ------> x4
 <<x^2>>   ------> x5
 <<y^2>>   ------> x6
  <<xy>>   ------> x7
 <<xp_x>>  ------> x8
 <<yp_y>>  ------> x9
 <<xp_y>>  ------> x10
 <<yp_x>>  ------> x11
<<p_xp_x>> ------> x12
<<p_yp_y>> ------> x13
<<p_xp_y>> ------> x14



Variable Matrix
------------------------------
x = [x1;  = x(1)
     x2;  = x(2)
     x3;  = x(3)
     x4;  = x(4)
     x5;  = x(5)
     x6;  = x(6)
     x7;  = x(7)
     x8;  = x(8)
     x9;  = x(9)
     x10;  = x(10)
     x11;  = x(11)
     x12;  = x(12)
     x13;  = x(13)
     x14]  = x(14)



Initial Conditions
-----------------------------
x0 = [0.625*factor;  = x1(0)
      0.325*factor;  = x2(0)
        0;           = x3(0)
        0;           = x4(0)
       0.5;          = x5(0)
       0.5;          = x6(0)
        0;           = x7(0)
        0;           = x8(0)
        0;           = x9(0)
        0;           = x10(0)
        0;           = x11(0)
       0.5;          = x12(0)
       0.5;          = x13(0)
        0];          = x14(0)
%}

%
tStart = tic;

% f and Initial Conditions
f  = @(t,x)([f1(t,x);f2(t,x);f3(t,x);f4(t,x);f5(t,x);f6(t,x);f7(t,x);f8(t,x);f9(t,x);f10(t,x);f11(t,x);f12(t,x);f13(t,x);f14(t,x)]);

f_classical = @(t2,x2)([fc1(t2,x2);fc2(t2,x2);fc3(t2,x2);fc4(t2,x2)]);

factor = 8;

x0 = [0.625*factor;
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
    
xc0 = [0.625*factor;
      0.325*factor;
        0;
        0];
    
% ode45
tspan = [0 3.75];
[t,x] = ode45(@appr, tspan, x0);
[t2,x2] = ode45(@appr2, tspan, xc0);


close all; % Closes the all figures
figure; 
%for f = 1: tspan = [0 30] and plot(t/15*3...
%for f = 2: tspan = [0 15] and plot(t/15*6...
%for f = 4: tspan = [0 7.5] and plot(t/15*12...
%for f = 8: tspan = [0 3.75] and plot(t/15*24...
plot(t/15*24,x(:,1),'LineWidth',2), grid; 
xlabel('{\boldmath$t/\tau_{L}$}'), ylabel('{\boldmath$<\hat{x}(t)>$}'); 
factorstr = sprintf('$f = %.0f$', factor);

title(factorstr);
%for f = 1: ylim([-0.8 0.8]);
%for f = 2: ylim([-1.5 1.5]);
%for f = 4: ylim([-2 3]);
%for f = 8: ylim([-3 6]);
ylim([-3 6]);
hold on

%for f = 1: tspan = [0 30] and plot(t/15*3...
%for f = 2: tspan = [0 15] and plot(t/15*6...
%for f = 4: tspan = [0 7.5] and plot(t/15*12...
%for f = 8: tspan = [0 3.75] and plot(t/15*24...
plot(t2/15*24,x2(:,1),'LineWidth',2), grid on; 

legend('Gaussian', 'Classic');

% set(gca,'LooseInset',get(gca,'TightInset'));
% saveas(gcf,'junk.eps')
% 
% set(findall(gcf,'Type','line'),'LineWidth',2)
% set(findall(gcf,'-property','FontSize'),'FontSize',12);

tEnd = toc(tStart)

function dxdt = appr(t,x)

dxdt = [x(3);
        x(4);
        -x(1)*(x(6)+x(2)^2)-2*x(7)*x(2);
        -x(2)*(x(5)+x(1)^2)-2*x(7)*x(1);
        2*x(8);
        2*x(9);
        x(11)+x(10);
        x(12)-x(5)*(x(6)+x(2)^2)-2*x(7)*(x(7)+x(1)*x(2));
        x(13)-x(6)*(x(5)+x(1)^2)-2*x(7)*(x(7)+x(1)*x(2));
        x(14)-x(7)*(x(5)+x(1)^2)-2*x(5)*(x(7)+x(1)*x(2));
        x(14)-x(7)*(x(6)+x(2)^2)-2*x(6)*(x(7)+x(1)*x(2));
        -2*x(8)*(x(6)+x(2)^2)-4*x(11)*(x(7)+x(1)*x(2));
        -2*x(9)*(x(5)+x(1)^2)-4*x(10)*(x(7)+x(1)*x(2));
        -x(10)*(x(6)+x(2)^2)-2*x(9)*(x(7)+x(1)*x(2))-x(11)*(x(5)+x(1)^2)-2*x(8)*(x(7)+x(1)*x(2))];
end
function dxdt2 = appr2(t2,x2)

dxdt2 = [x2(3);
        x2(4);
        -x2(1)*(x2(2)^2);
        -x2(2)*(x2(1)^2)];
end