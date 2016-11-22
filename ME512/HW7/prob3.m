clear; clc; close all;

R = 0.3; % m
omega = 2*pi; % rad/s
nu = 1.004*10^-6;

r = 0 : 0.01 : R;
%t = 0 : 500 : 50000;
%t = 450;
t = [450,3600,1000000];

that = t*nu/R^2;

rhat = r/R;

n_terms = 20;

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
            
            term(n+1) = 2*besselj(1,l_n(n+1)*rhat(j))/(l_n(n+1)*besselj(0,l_n(n+1)))...
                *exp(-l_n(n+1)^2.0*that(i));
            
        end
        
        v(j,i) = rhat(j) + sum(term); 
        
    end
    
end

v_theta = R*omega*v;

hold on
plot(r,v(:,1),'--');
plot(r,v(:,2),'--x');
plot(r,v(:,3));
legend('t = 450s','t = 3600s','steady-state');
xlabel('radial position (m)');
ylabel('V_{normalized} (m/s)');

