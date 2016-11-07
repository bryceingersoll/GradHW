clear; clc; close all;

rho = 999; % kg/m^3
g = 9.81; % m/s^2

%choose omega, gamma arbitrarily
w = 1;
gamma = 1;

p_atm = 1;

r = 0.2 : 0.01 : 1.0;

h_force = zeros(length(r),1);
h_free = zeros(length(r),1);

for i = 1 : length(r)
   
    h_force(i) = ((w*r(i))^2)/(2*g);
    
    h_free(i) = -gamma^2/(2*pi^2*g*r(i)^2);
    
end

hold on
plot(r,h_force,'--')
plot(r,h_free)
legend('Forced Vortex Height','Free Vortex Height');
title('Height Comparison of Free Surface');
xlabel('Radial Position (m)');
ylabel('Free Surface Shape (m)');