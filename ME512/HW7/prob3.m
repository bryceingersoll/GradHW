clear; clc; close all;

R = 0.3; % m
omega = 2*pi; % rad/s
nu = 1*10^-6;

r = 0 : 0.001 : R;
%t = 0 : 500 : 50000;
%t = 360;
t= 4500;

n_terms = 10;

lambda = 0.0001 : 0.0001 : 100;

f = zeros(length(lambda),1);

for i = 1 : length(lambda)
    
    f(i) = besselj(1,lambda(i));
end

l_n = zeros(n_terms,1);
ng = 1;
for i = 2 : length(lambda)
   
    if f(i) > 0 && f(i-1) < 0 % neg to pos
        l_n(ng) = lambda(i);
        ng = ng + 1;
    end
       if f(i) < 0 && f(i-1) > 0 % pos to neg
        l_n(ng) = lambda(i);
        ng = ng + 1;
       end
end

%l_n = [0; l_n(1:n_terms-1)];

v = zeros(length(r),length(t));

term = zeros(n_terms,1);

for i = 1 : length(t)
    
    for j = 1 : length(r)
        
        for n = 0 : n_terms-1
            
            term(n+1) = 2*besselj(1,l_n(n+1)*r(j)/R)/(l_n(n+1)*besselj(0,l_n(n+1)))...
                *exp(-l_n(n+1)^2.0*t(i)*nu/R^2);
            
        end
        
        v(j,i) = r(j) + sum(term); 
        
    end
    
end

v_theta = R*omega*v;

plot(r,v(:,1));