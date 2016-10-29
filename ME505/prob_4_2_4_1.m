clear; clc; close all;

% boundaries
% M = 1; % x boundary
% L = 1; % y boundary

x = linspace(0,1,20);
y = linspace(0,1,20);

% number of terms
n = 10;

% function
f = zeros(length(y),1);
for i = 1 : length(y)
    
    f(i) = y(i)*(1-y(i));
    
end

%terms
u = zeros(length(x),length(y));
un = zeros(length(x),length(y), n);
Xn = zeros(length(x),1);
Yn = zeros(length(y),1);
Cn = zeros(length(y),1);
top = Cn;
bot = Cn;

for i = 1 : n
    
    for j = 1 : length(x)
        
        %define Xn
        Xn(j) = cosh((n+1)*pi*(x(j)-1));
        
        for k = 1 : length(y)
            
            %define Yn
            Yn(k) = sin(pi*y(k)*(n+1));
            
            %define Cn
        
            fun = @(y) y.*(1-y).*sin((n+1)*pi.*y)*2;
            
            top(k) = integral(fun, 0, 1);
            
            fun = @(y) (sin(pi.*y.*(n+1))).^2;
            
            bot(k) = cosh(pi*(n+1))*integral(fun,0,1);
            
            Cn(k) = top(k)/bot(k);
            
            % calculate u
            un(j,k,n) = Xn(j)*Yn(k)*Cn(k);
            
        end
        
    end
    
    % add terms
    u = u + un(:,:,n);
    
end

surf(x,y,u);
