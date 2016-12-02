clear; clc; close all;

% not sure
a = 1;
ui = 1;
nu = 1;

%fun = @(x) exp(-x.^2).*log(x).^2;
%q = integral(fun,0,Inf);

x = 0.1 : 0.001 : pi/2.0;

fun = @(x) (ui*(1.814.*x/a-0.271*(x/a).^3.0-0.0471*(x/a).^5.0)).^5.0;

theta_squared = zeros(length(x),1);
theta = zeros(length(x),1);
lambda = zeros(length(x),1);

for i = 1 : length(x)
   
    theta_squared(i) = nu*0.45/((ui*(1.814.*x(i)/a-0.271*(x(i)/a).^3.0-0.0471*(x(i)/a)^5.0)))^6.0...
        *integral(fun,0,x(i));
    
    theta(i) = theta_squared(i)^0.5;
    
    lambda(i) = theta_squared(i)/nu*ui*(1.814/a-3*0.271*x(i).^2.0/a.^3.0-5*0.0471*x(i).^4.0/a.^5.0);
    
end

plot(x,lambda);