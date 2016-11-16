clear; clc; close all;

U = 0.3; % m/s
h = 0.02; % m
nu = 1.983*10^-5/1.225; % m^2/s

%y = 0 : 0.0001 : h;
y = h/2.0;

t = 7.3 : 0.0001 : 7.4;


n_terms = 20;

u1 = zeros(length(y),length(t));

u1n = zeros(length(y),length(t),n_terms);

u = zeros(length(y),length(t));

for i = 1 : length(t)
    
    for j = 1 : length(y)
        for n = 1 : n_terms
            
            u1n(j,i,n) = -U*2/pi/n*sin(n*pi*y(j)/h)*exp(-n^2*pi^2*t(i)*nu/h^2);
            
        end
        u1(j,i) = sum(u1n(j,i,:));
    end
       
    u(:,i) = U*y'/h + u1(:,i);
end

hold on;
plot(u(:,10),y);
plot(u(:,20),y);
plot(u(:,30),y);
plot(u(:,40),y);
