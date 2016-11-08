clear; clc; close all;

L = 2;
M = 1;

n_terms = 5;


x = 0 : L/10.0 : L;
y = 0 : M/10.0 : M;

%Poisson Equation
X_m = zeros(n_terms,1);
Y_n = zeros(n_terms,1);
mu_m = zeros(n_terms,1);
rho_n = zeros(n_terms,1);

A_mn = zeros(n_terms);

for i = 1 : length(x)
    
    for j = 1 : length(y)
        
        for m = 1 : n_terms
            
            %Define X_m
            X_m(m) = sin((m+1)*pi/L*x(i));
            
            %Define mu_m
            mu_m(m) = ((m+1)*pi/L)^2.0;
            
            for n = 1 : n_terms
                
                %Define Y_n
                Y_n(n) = sin((2*n+1)*pi/(2*M)*y(j));
                
                %Define rho_n
                rho_n(n) = ((2*n+1)*pi/(2*M))^2.0;
                
                %Define A_mn
                fun = @(x,y) 0.1 .* x .* y .* (L - x) .* (M - y) .* sin((m+1).*pi./L.*x) .*sin((2.*n+1)*pi./(2.*M).*y);
                
                top = integral2(fun,0,L,0,M);
                
                fun_x = @(x)  sin((m+1).*pi./L.*x).^2.0;
                
                fun_y = @(y) sin((2.*n+1)*pi./(2.*M).*y).^2.0;
                
                int_x = integral(fun_x,0,L);
                
                int_y = integral(fun_y,0,M);
                
                bot = int_x*int_y*(mu_m(m) + rho_n(n));
                
                A_mn(m,n) = top/bot;
                
                
            end
            
        end
        
    end
    
end