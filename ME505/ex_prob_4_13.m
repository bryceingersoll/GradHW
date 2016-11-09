clear; clc; close all;

tic

%L = 2;
%M = 1;
L = 6;
M = 4;

n_terms = 10;


x = 0 : L/5.0 : L;
y = 0 : M/5.0 : M;

%Poisson Equation
X_m = zeros(n_terms,1);
Y_n = zeros(n_terms,1);
mu_m = zeros(n_terms,1);
rho_n = zeros(n_terms,1);
A_mn = zeros(n_terms);
v_mn = zeros(length(x),length(y),n_terms, n_terms);
v = zeros(length(x),length(y));

for i = 1 : length(x)
    
    for j = 1 : length(y)
        
        for m = 0 : n_terms-1
            
            %Define X_m
            %X_m(m+1) = sin((m+1)*pi/L*x(i));
            X_m(m+1) = cos((2*m+1)*pi/(2*L)*x(i));
            
            %Define mu_m
            %mu_m(m+1) = ((m+1)*pi/L)^2.0;
            mu_m(m+1) = ((2*m+1)*pi/(2*L))^2.0;
            
            for n = 0 : n_terms-1
                
                %Define Y_n
                Y_n(n+1) = sin((2*n+1)*pi/(2*M)*y(j));
                
                %Define rho_n
                rho_n(n+1) = ((2*n+1)*pi/(2*M))^2.0;
                
                %Define A_mn
                %fun = @(x,y) 0.1 .* x .* y .* (L - x) .* (M - y) .* sin((m+1).*pi./L.*x) .*sin((2.*n+1)*pi./(2.*M).*y);
                fun = @(x,y) 0.02 .* x .* y .* (L - x) .* (M - y) .* sin((2.*m+1).*pi./(2.*L).*x) .*sin((2.*n+1)*pi./(2.*M).*y);
                
                top = integral2(fun,0,L,0,M);
                
                fun_x = @(x)  sin((2.*m+1).*pi./(2.*L).*x).^2.0;
                
                fun_y = @(y) sin((2.*n+1)*pi./(2.*M).*y).^2.0;
                
                int_x = integral(fun_x,0,L);
                
                int_y = integral(fun_y,0,M);
                
                bot = int_x*int_y*(mu_m(m+1) + rho_n(n+1));
                
                A_mn(m+1,n+1) = -top/bot;
                %A_mn(m+1,n+1) = top*(-16*L*M)/(4*M^2*(m+1)^2+L^2*(2*n+1)^2)/pi^2;
                
                %calculate v_mn
                v_mn(i,j,m+1,n+1) = X_m(m+1)*Y_n(n+1)*A_mn(m+1,n+1);
                
            end
            
        end
        
    end
    
end

v = zeros(length(x),length(y));

for m = 1 : n_terms
    
    for n = 1 : n_terms
        
        v = v + v_mn(:,:,m,n);
        
    end
    
end

toc

% % Laplace Equation
% %terms
% w = zeros(length(x),length(y));
% wn = zeros(length(x),length(y), n_terms);
% Xn = zeros(length(x),1);
% Yn = zeros(length(y),1);
% Cn = zeros(length(y),1);
% top = Cn;
% bot = Cn;
% 
% for n = 0 : n_terms-1
%     
%     for j = 1 : length(x)
%         
%         %define Xn
%         Xn(j) = cosh((n+1)*pi*(x(j)-1));
%         
%         for k = 1 : length(y)
%             
%             %define Yn
%             Yn(k) = sin(pi*y(k)*(n+1));
%             
%             %define Cn
%         
%             fun = @(y) y.*(1-y).*sin((n+1)*pi.*y)*2;
%             
%             top(k) = integral(fun, 0, 1);
%             
%             fun = @(y) (sin(pi.*y.*(n+1))).^2;
%             
%             bot(k) = cosh(pi*(n+1))*integral(fun,0,1);
%             
%             Cn(k) = top(k)/bot(k);
%             
%             % calculate u
%             un(j,k,n+1) = Xn(j)*Yn(k)*Cn(k);
%             
%         end
%         
%     end
%     
%     % add terms
%     u = u + un(:,:,n+1);
%     
% end

surf(x,y,v);
xlabel('x');
ylabel('y');

