clear; clc; close all;
tic
% boundaries
% M = 1; % x boundary
% L = 1; % y boundary

x = linspace(0,1,20);
y = linspace(0,1,20);

% number of terms
n_terms = 10;

%terms
u = zeros(length(x),length(y));
un = zeros(length(x),length(y), n_terms);
Xn = zeros(length(x),1);
Yn = zeros(length(y),1);
Cn = zeros(length(y),1);
top = Cn;
bot = Cn;

lambda = .001 : 0.00001 : 500;

function_value = zeros(length(lambda),1);
for i = 1 : length(lambda)
    
    function_value(i) = sin(lambda(i)) + lambda(i)*cos(lambda(i));
    
end

lambda_guess = zeros(n_terms,1);
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

lambda_guess = lambda_guess(1:n_terms,1);

%could use FZERO, but this is good enough
eigen_values = lambda_guess;

for n = 0 : n_terms-1
    
    for j = 1 : length(x)
        
        %define Xn
        Xn(j) = sin(eigen_values(n+1)*x(j));
        
        for k = 1 : length(y)
            
            %define Yn
            Yn(k) = sinh(eigen_values(n+1)*(y(k)-1));
            
            %define Cn
        
            fun = @(x) x.*(1.0-2.0*x/3.0).*sin(x.*eigen_values(n+1));
            
            top(k) = integral(fun, 0, 1);
            
            fun = @(x) sin(x.*eigen_values(n+1)).^2;
            
            bot(k) = sinh(-eigen_values(n+1))*integral(fun,0,1);
            
            Cn(k) = top(k)/bot(k);
            
            % calculate u
            un(j,k,n+1) = Xn(j)*Yn(k)*Cn(k);
            
        end
        
    end
    
    % add terms
    u = u + un(:,:,n+1);
    
end
toc
surf(x,y,u);
xlabel('y');
ylabel('x');