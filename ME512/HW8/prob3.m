clear; clc; close all;

tic

% not sure
a = 1;
ui = 1;
nu = 1;
mu = 1;
rho = 1;

check = 0;

%fun = @(x) exp(-x.^2).*log(x).^2;
%q = integral(fun,0,Inf);

x = 0.1 : 0.0001 : pi/2.0;

fun = @(x) (ui*(1.814.*x/a-0.271*(x/a).^3.0-0.0471*(x/a).^5.0)).^5.0;

theta_squared = zeros(length(x),1);
theta = zeros(length(x),1);
lambda = zeros(length(x),1);
S = zeros(length(x),1);
tau_w = zeros(length(x),1);
c_f = zeros(length(x),1);

for i = 1 : length(x)
   
    %calculate theta squared
    theta_squared(i) = nu*0.45/((ui*(1.814.*x(i)/a-0.271*(x(i)/a).^3.0-0.0471*(x(i)/a)^5.0)))^6.0...
        *integral(fun,0,x(i));
    
    %calculate theta
    theta(i) = theta_squared(i)^0.5;
    
    %calculate lambda
    lambda(i) = theta_squared(i)/nu*ui*(1.814/a-3*0.271*x(i).^2.0/a.^3.0-5*0.0471*x(i).^4.0/a.^5.0);
    
    if lambda(i) < -0.09 && check == 0
        x_sep = x(i);
        
        x_sep*180/pi
        
        check = 1;
    end
        
    %calculate S(lambda)
    S(i) = (lambda(i)+0.9)^(0.62);
    
    %calculate tau_w
    tau_w(i) = mu*(ui*(1.814.*x(i)/a-0.271*(x(i)/a).^3.0-0.0471*(x(i)/a).^5.0))/theta(i)*S(i);
    
    %calculate C_f
    c_f(i) = 2*tau_w(i)/(rho*(ui*(1.814.*x(i)/a-0.271*(x(i)/a).^3.0-0.0471*(x(i)/a).^5.0))^2.0);
    
end

figure(1);
hold on;
xlabel('Distance along outside of cylinder');
ylabel('Wall Shear Stress');
title('\tau_w vs. Distance');
plot(x,tau_w);
hold off;

figure(2);
hold on;
xlabel('Distance along outside of cylinder');
ylabel('Boundary layer thickness');
title('\theta vs. Distance');
plot(x,theta);
hold off;

toc;