clear; clc; close all;

g = 9.81; % m/s
gx = g*sind(20.0);

h = 1;

mu = 8.9*10^-4;
rho = 999;
nu = mu/rho;

umax = 1.0;

%t = 0 : 0.01 : 1.0;
%t = 1000;

y = 0 : 0.01 : h;

yhat = y/h;
ghat = gx*h^2/nu/umax;
%that = t*nu/h^2;
that = 0.2365;
t = that*h^2/nu;
n_terms = 10;

uhat = zeros(length(y),1);

%U = g/2(y-y^2)+SUM2g/(n.pi)^3 (cos(n.pi)-1)sing(n.pi.y)exp(-n.pi)^2.t
term = zeros(length(y),n_terms);

for i = 1 : length(y)
    for n = 1 : n_terms
        
        term(i,n) = 2*ghat/(n*pi)^3*(cos(n*pi)-1)*sin(n*pi*yhat(i))*exp(-n^2*pi^2*that);
        
    end
    
    uhat(i) = ghat*(yhat(i)-yhat(i)^2) + sum(term(i,:));
    
end

plot(uhat(1:51)/ghat,y(1:51));
    