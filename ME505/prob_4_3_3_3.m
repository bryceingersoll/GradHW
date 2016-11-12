clear; clc; close all;

tic

L = 2;
M = 1;

n_terms = 10;


x = 0 : L/20.0 : L;
y = 0 : M/10.0 : M;

pe = 1;
le = 1;

if pe == 1
    %Poisson Equation
    X_m = zeros(n_terms,1);
    Y_n = zeros(n_terms,1);
    mu_m = zeros(n_terms,1);
    rho_n = zeros(n_terms,1);
    A_mn = zeros(n_terms);
    v_mn = zeros(length(x),length(y),n_terms, n_terms);
    v = zeros(length(x),length(y));
    
    for m = 0 : n_terms-1
        
        %Define mu_m
        mu_m(m+1) = ((m+1)*pi/L)^2.0;
        
        for n = 0 : n_terms-1
            
            %Define rho_n
            rho_n(n+1) = (pi*n/M)^2.0;
            
            %Define A_mn
            fun = @(x,y) 0.1 .* x .* y .* (L - x) .* (M - y) .* sin((m+1).*pi./L.*x)...
                .*cos(pi*n.*y);
            
            top = integral2(fun,0,L,0,M);
            
            fun_x = @(x)  sin((m+1).*pi./(L).*x).^2.0;
            
            fun_y = @(y) cos(pi*n/M.*y).^2.0;
            
            int_x = integral(fun_x,0,L);
            
            int_y = integral(fun_y,0,M);
            
            bot = int_x*int_y*(mu_m(m+1) + rho_n(n+1));
            
            A_mn(m+1,n+1) = -top/bot;
            %A_mn(m+1,n+1) = top*(-16*L*M)/(4*M^2*(m+1)^2+L^2*(2*n+1)^2)/pi^2;
            
            
            for i = 1 : length(x)
                
                %Define X_m
                X_m(m+1) = sin((m+1)*pi/L*x(i));
                
                
                
                for j = 1 : length(y)
                    
                    %Define Y_n
                    Y_n(n+1) = cos(pi*n/M*y(j));
                    
                    %calculate v_mn
                    v_mn(i,j,m+1,n+1) = X_m(m+1)*Y_n(n+1)*A_mn(m+1,n+1);
                    
                end
                
            end
            
            v = v + v_mn(:,:,m+1,n+1);
            
        end
        
    end
    
    v = zeros(length(x),length(y));
    
    for m = 1 : n_terms
        
        for n = 1 : n_terms
            
            v = v + v_mn(:,:,m,n);
            
        end
        
    end
    
end

if le == 1
    % 1st Laplace Equation
    %terms
    w1 = zeros(length(x),length(y));
    wn = zeros(length(x),length(y), n_terms);
    Xn = zeros(length(x),1);
    Yn = zeros(length(y),1);
    Bn = 0;
    top = Bn;
    bot = Bn;
    
    for n = 0 : n_terms-1
        
        %define Bn
        
        if n ==0
            fun = @(y) y;
            
            Bn = -integral(fun,0,M)/(M*L);
            
        else
        fun = @(y) y.*cos(n*pi.*y./M);
        
        top = integral(fun, 0, M);
        
        %fun = @(y) (sin(pi.*y.*(2*n+1)/(2*M))).^2;
        
        bot = M*sinh(n*pi/M*L)*0.5;
        
        Bn = -top/bot;
        end
        
        for j = 1 : length(x)
            
            %define Xn
            if n == 0
                Xn(j) = x(j)-L;
            else
            Xn(j) = sinh(n*pi/M*(x(j)-L));
            end
            
            for k = 1 : length(y)
                
                %define Yn
                Yn(k) = cos(n*pi/M*y(k));
                
                % calculate u
                wn(j,k,n+1) = Xn(j)*Yn(k)*Bn;
                
            end
            
        end
        
        % add terms
        w1 = w1 + wn(:,:,n+1);
        
    end
    
        % Laplace Equation
    %terms
    w2 = zeros(length(x),length(y));
    wn = zeros(length(x),length(y), n_terms);
    Xn = zeros(length(x),1);
    Yn = zeros(length(y),1);
    Dn = 0;
    top = Dn;
    bot = Dn;
    
    for n = 0 : n_terms-1
        
        %define Bn
        
        if n == 0
            fun = @(y) y.*(M-y);
            
            Dn = integral(fun,0,M)/(L*M);
            
        else
            
        fun = @(y) y.*(M-y).*cos(n*pi/M.*y);
        
        top = integral(fun, 0, M);
        
        fun = @(y) cos(n*pi/M.*y).^2;
        
        bot = sinh(n*pi/M*L)*integral(fun,0,1);
        
        Dn = top/bot;
        end
        
        for j = 1 : length(x)
            
            %define Xn
            if n == 0
                Xn(j) = x(j);
            else
            Xn(j) = sinh(n*pi/M*x(j));
            end
            
            for k = 1 : length(y)
                
                %define Yn
                Yn(k) = cos(n*pi/M*y(k));
                
                % calculate u
                wn(j,k,n+1) = Xn(j)*Yn(k)*Dn;
                
            end
            
        end
        
        % add terms
        w2 = w2 + wn(:,:,n+1);
        
    end
    
    
end

toc

if le == 1 && pe == 1
    
    u = v + w1 + w2;
elseif le == 0
    u = v;
else
    u = w1 + w2;
end

surf(x,y,u');
xlabel('x');
ylabel('y');


