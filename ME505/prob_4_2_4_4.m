clear; clc; close all;

% boundaries
% r = 1; % radial boundary
% theta = pi/3; % theta boundary

r = linspace(0,1,20);
theta = linspace(0,pi/3,20);

% number of terms
n_terms = 10;

%terms
u = zeros(length(r),length(theta));
un = zeros(length(r),length(theta), n_terms);
Rn = zeros(length(r),1);
Thetan = zeros(length(theta),1);
Cn = zeros(length(theta),1);
top = Cn;
bot = Cn;

for n = 0 : n_terms-1
    
    for j = 1 : length(r)
        
        %define Rn
        Rn(j) = r(j)^(1.5*(2*n+1));
        
        for k = 1 : length(theta)
            
            %define Thetan
            Thetan(k) = cos(1.5*(2*n+1)*theta(k));
            
            %define Cn
        
            fun = @(theta) (1-9/pi^2.*theta.^2).*cos(1.5*(2.*n+1).*theta);
            
            top(k) = integral(fun, 0, pi/3);
            
            fun = @(theta) cos(1.5*(2.*n+1).*theta).^2;
            
            bot(k) = integral(fun,0,pi/3);
            
            Cn(k) = top(k)/bot(k);
            
            % calculate u
            un(j,k,n+1) = Rn(j)*Thetan(k)*Cn(k);
            
        end
        
    end
    
    % add terms
    u = u + un(:,:,n+1);
    
end

%transform from polar to cartesian

[R,Theta] = meshgrid(r,theta);

surf(R.*cos(Theta),R.*sin(Theta),u');
