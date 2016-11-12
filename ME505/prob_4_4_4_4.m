clear; clc; close all;

tic

L = 1;
M = pi;

n_terms = 10;

k = 0.1;

r = 0 : L/20.0 : L;
theta = 0 : pi/20.0 : pi;

%set time
t = 0.1;

Theta_m = zeros(n_terms,1);
R_mn = zeros(n_terms);
rho_m = zeros(n_terms,1);
mu_mn = zeros(n_terms);
T_mn = zeros(n_terms);
C_mn = zeros(n_terms);
u_mn = zeros(length(r),length(theta),n_terms, n_terms);
u = zeros(length(r),length(theta));

for m = 0 : n_terms-1
    
    %Define rho_m
    rho_m(m+1) = m+1;
    
    %Define mu_mn
    p_mu = 0.001 : 0.001 : 100;
    
    bessel_value = besselj(m+1,p_mu);
    
    ng = 1;
    
    for i = 2 : length(bessel_value)
        
        if bessel_value(i) < 0 && bessel_value(i-1) > 0 %goes pos to neg
            if ng < n_terms+1
            mu_mn(m+1,ng) = p_mu(i);
            ng = ng + 1;
            end
        end
        
        if bessel_value(i) > 0 && bessel_value(i-1) < 0 %goes neg to pos
            
            if ng < n_terms+1
            mu_mn(m+1,ng) = p_mu(i);
            ng = ng + 1;
            end
        end
        
    end
    
    
    for n = 0 : n_terms-1
        
        
        %mu_mn(n+1,m+1) = (n+1)*pi/2;
        
        %Define C_mn
        fun = @(r,theta) (r-r.^3).*sin((m+1).*theta).*besselj(m+1,mu_mn(m+1,n+1).*r).*sin((m+1).*theta).*r;
        
        top = integral2(fun,0,1,0,pi);
        
        fun_r = @(r)  besselj(m+1,mu_mn(m+1,n+1).*r).^(2.0).*r;
        
        fun_theta = @(theta) sin((m+1).*theta).^2.0;
        
        int_r = integral(fun_r,0,L);
        
        int_theta = integral(fun_theta,0,pi);
        
        bot = int_theta*int_r;
        
        C_mn(m+1,n+1) = top/bot;
        
        %Define Tmn
        T_mn(m+1,n+1) = exp(-mu_mn(m+1,n+1).^2.0*pi*k*t);
        
        for i = 1 : length(r)
            
            %Define R_mn
            R_mn(m+1,n+1) = besselj(m+1,mu_mn(m+1,n+1)*r(i));
            
            for j = 1 : length(theta)
                
                %Define Theta_m
                Theta_m(n+1) = sin((m+1)*theta(j));
             
                %calculate u_mn
                u_mn(i,j,m+1,n+1) = T_mn(m+1,n+1)*Theta_m(m+1)*R_mn(m+1,n+1)*C_mn(m+1,n+1);
                
            end
            
        end
        
        u = u + u_mn(:,:,m+1,n+1);
        
    end
    
end

u = zeros(length(r),length(theta));

for m = 1 : n_terms
    
    for n = 1 : n_terms
        
        u = u + u_mn(:,:,m,n);
        
    end
    
end

toc

[R,Theta] = meshgrid(r,theta);

surf(R.*cos(Theta),R.*sin(Theta),u');
