clear; clc; close all;

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
Xn_a = zeros(length(x),1);
Yn_a = zeros(length(y),1);
An = zeros(length(y),1);

Xn_b = zeros(length(x),1);
Yn_b = zeros(length(y),1);
Bn = zeros(length(y),1);

top = An;
bot = An;

for n = 0 : n_terms-1
    
    for j = 1 : length(x)
        
        %define Xn
        Xn_a(j) = sin((2*n+1)*pi*x(j)/2);
        
        Xn_b(j) = cosh((n+1)*pi*(x(j)-1));
        
        for k = 1 : length(y)
            
            %define Yn
            Yn_a(k) = sinh((2*n+1)*pi*y(k)/2);
            
            Yn_b(k) = sin((n+1)*pi*y(k));
            
            %define An
        
            fun = @(x) x.*(1-0.5*x).*sin((2*n+1)*pi.*x./2);
            
            top(k) = 2*integral(fun, 0, 1);
            
            bot(k) = sinh((2*n+1)*pi/2);
            
            An(k) = top(k)/bot(k);
            
            %define Bn
            fun = @(y) y.*(1-y).*sin((n+1)*pi.*y);
            
            top(k) = 2*integral(fun, 0, 1);
            
            bot(k) = cosh((n+1)*pi);
            
            Bn(k) = top(k)/bot(k);
            
            % calculate u
            un(j,k,n+1) =  Xn_a(j)*Yn_a(k)*An(k) +Xn_b(j)*Yn_b(k)*Bn(k);
            
        end
        
    end
    
    % add terms
    u = u + un(:,:,n+1);
    
end

surf(x,y,u);
xlabel('y');
ylabel('x');
