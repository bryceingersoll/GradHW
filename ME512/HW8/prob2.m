clear; clc; close all;

tic

% for water, arbitrary

mu = 1.002*10^-3;
rho = 999.0;

H = 1;

% example
% syms x
% eqn = sin(x) == 1;
% solx = solve(eqn,x)

u0 = 1.0;

scale = 2500;

x = scale*(0.0 : 0.2 : 10.0);

sol_uc = zeros(length(x),1);
delta = zeros(length(x),1);

for i = 1 : length(x)
    
    syms uc
    eqn = 3*H - 3*u0*H/uc == (60*x(i)*mu/rho/uc)^0.5;
    sol_uc(i) = solve(eqn,uc);
    
    
    delta(i) = (60*x(i)*mu/rho/sol_uc(i))^0.5;
    
    y = 0 : delta(i)/100.0 : delta(i);
    
    ub = zeros(length(y),1);
    
    for j = 1 : length(y)
    
        ub(j) = sol_uc(i)*(2*y(j)/delta(i) - (y(j)/delta(i))^2.0);
    
    end
    
    figure(1);
    hold on
    plot(ub,y,'b');
    plot(ub,2*H*ones(length(y),1)'-y,'b');
    plot([sol_uc(i), sol_uc(i)], [delta(i), 2*H-delta(i)],'b');
    xlabel('u/u_0');
    ylabel('y (H)');
    xlim([0,1.6]);
    ylim([0,2]);
    a = [0.6 0.9];
    b = [0.5 0.5];
    annotation('textarrow',a,b,'String',' Fluid Velocity Profile as x increases')
    hold off;
    
    
    
    
end

l_e = H^2*rho*u0/(40*mu);

figure(2);
hold on;
plot(x/scale,delta);
xlabel('distance from start of tube');
ylabel('\delta/H')
title('\delta(x) for flow between parallel plates');
ylim([0,1]);
hold off;
% 
% figure(2);
% plot(x/scale,sol_uc);

toc