clear; clc; close all;
tic
% boundaries
% M = 1; % x boundary
% L = 1; % y boundary

x = linspace(0,1,40);
y = linspace(0,1,40);

% number of terms
n_terms = 10;

%terms
u = zeros(length(x),length(y));
un = zeros(length(x),length(y), n_terms);
Xn = zeros(length(x),1);
Yn = zeros(length(y),1);
Cn = zeros(length(y),1);
top = Cn;
bot = Cn;

eigen_values = solve_eigen_values(n_terms);


for n = 0 : n_terms-1
    
    for j = 1 : length(x)
        
        %define Xn
        Xn(j) = sin(eigen_values(n+1)*x(j));
        
        for k = 1 : length(y)
            
            %define Yn
            Yn(k) = sinh(eigen_values(n+1)*(y(k)-1));
            
            %define Cn
        
            fun = @(x) x.*(1.0-2.0*x/3.0).*sin(x.*eigen_values(n+1));
            
            top(k) = integral(fun, 0, 1);
            
            fun = @(x) sin(x.*eigen_values(n+1)).^2;
            
            bot(k) = sinh(-eigen_values(n+1))*integral(fun,0,1);
            
            Cn(k) = top(k)/bot(k);
            
            % calculate u
            un(j,k,n+1) = Xn(j)*Yn(k)*Cn(k);
            
        end
        
    end
    
    % add terms
    u = u + un(:,:,n+1);
    
end
toc
surf(x,y,u);
xlabel('x');
ylabel('y');