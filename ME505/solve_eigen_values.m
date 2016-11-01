%clear; clc; close all;
function eigen_values = solve_eigen_values(n)
%n = 10;
% for problem 4.2.4:2
lambda = .001 : 0.00001 : 500;

function_value = zeros(length(lambda),1);
for i = 1 : length(lambda)
    
    function_value(i) = lambda(i)*sin(lambda(i))+cos(lambda(i));
    
end

lambda_guess = zeros(n,1);
ng = 1;

for i = 2 : length(function_value)
    
    if function_value(i) < 0 && function_value(i-1) > 0 %goes pos to neg
        lambda_guess(ng) = lambda(i);
        ng = ng + 1;
    end
    
    if function_value(i) > 0 && function_value(i-1) < 0 %goes neg to pos
        lambda_guess(ng) = lambda(i);
        ng = ng + 1;
    end
    
end

lambda_guess = lambda_guess(1:n,1);
eigen_values = lambda_guess;

% hold on
% plot(lambda,function_value)
% plot(lambda_guess,zeros(n,1),'*')

% lambda = linspace(0,100,1000);
%
% hold on
% plot(lambda,lambda.*sin(lambda)+cos(lambda))
% plot(lambda,zeros(length(lambda),1))
