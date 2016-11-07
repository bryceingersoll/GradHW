clear; clc; close all;

% boundaries
% r = 1; % radial boundary
% theta = pi/3; % theta boundary

r = linspace(1,2,20);
theta = linspace(0,pi/2,20);

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
        Rn(j) = sin((2*n+1)*pi/(2*log(2))*log(r(j)));
        
        for k = 1 : length(theta)
            
            %define Thetan
            Thetan(k) = sinh((2*n+1)*pi/(2*log(2))*(theta(k)-pi/2));
            
            %define Cn
        
            fun = @(r) (-r.^2+4.*r-3).*sin((2.*n+1).*pi./(2.*log(2)).*log(r))./r;
            
            top(k) = integral(fun,1,2);
            
            bot(k) = -sinh((2*n+1)*pi^2/(4*log(2)))*log(2)*0.5;
            
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
