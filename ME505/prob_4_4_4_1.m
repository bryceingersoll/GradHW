clear; clc; close all;

tic

L = 1;
M = 1;

n_terms = 10;

k = 0.1;

x = 0 : L/40.0 : L;
y = 0 : M/40.0 : M;

%set time
t = 0.01;

X_m = zeros(n_terms,1);
Y_n = zeros(n_terms,1);
alpha_m = zeros(n_terms,1);
beta_n = zeros(n_terms,1);
T_mn = zeros(n_terms);
C_mn = zeros(n_terms);
u_mn = zeros(length(x),length(y),n_terms, n_terms);
u = zeros(length(x),length(y));

for m = 0 : n_terms-1
    
    %Define mu_m
    alpha_m(m+1) = pi*(m+1);
    
    for n = 0 : n_terms-1
        
        %Define rho_n
        beta_n(n+1) = pi*(n+1);
        
        %Define A_mn
        fun = @(x,y) x.*(1-x).*(1-y).*sin((m+1)*pi.*x).*sin((n+1)*pi.*y);
        
        top = integral2(fun,0,L,0,M);
        
        fun_x = @(x)  sin((m+1)*pi.*x).^2.0;
        
        fun_y = @(y) sin((n+1)*pi.*y).^2.0;
        
        int_x = integral(fun_x,0,L);
        
        int_y = integral(fun_y,0,M);
        
        bot = int_x*int_y;
        
        %C_mn(m+1,n+1) = top/bot;
        C_mn(m+1,n+1) = top*4.0;
        %A_mn(m+1,n+1) = top*(-16*L*M)/(4*M^2*(m+1)^2+L^2*(2*n+1)^2)/pi^2;
        
        %Define Tmn
        T_mn(m+1,n+1) = exp(-(alpha_m(m+1)^2.0 + beta_n(n+1)^2.0)*pi*k*t);
        
        for i = 1 : length(x)
            
            %Define X_m
            X_m(m+1) = sin(pi*x(i)*(m+1));
            
            
            
            for j = 1 : length(y)
                
                %Define Y_n
                Y_n(n+1) = sin(pi*y(j)*(n+1));
                
                
                
                %calculate u_mn
                u_mn(i,j,m+1,n+1) = T_mn(m+1,n+1)*X_m(m+1)*Y_n(n+1)*C_mn(m+1,n+1);
                
            end
            
        end
        
        u = u + u_mn(:,:,m+1,n+1);
        
    end
    
end

u = zeros(length(x),length(y));

for m = 1 : n_terms
    
    for n = 1 : n_terms
        
        u = u + u_mn(:,:,m,n);
        
    end
    
end

toc

surf(x,y,u');
xlabel('x');
ylabel('y');