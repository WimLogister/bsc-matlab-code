function [ dx ] = dosedyn( t,x,dosage,p )
% Population and strategy dynamics for a single scalar phenotype strategy
% according to Brown et al (2015).
% x is the population, dosage is the treatment intensity, p is a struct
% array storing system parameters

    % Stores dx/dt and du/dt
    dx=zeros(2,1);
    
    % Series of checks to keep population size and resistance amount within
    % reasonable bounds
    if x(1) > 100
        x(1) = 100;
    end
    
    if x(1) < 0
        x(1) = 0;
    end
    
    if x(2) < 0
        x(2) = 0;
    end
    
    % Compute carrying capacity K
    K=p.Kmax.*exp((-x(2).^2)./(2*p.sig^2));
    
    % K==0 leads to numerical problems (NaN)
    if K == 0
        K = eps;
    end
    
    % Compute numerator of treatment efficacy mu
    top = (dosage*p.m*p.N^p.alpha)/p.N;
    
    % Compute denominator of mu
    % Note: since all models we consider only use a single scalar strategy, u = v
    % and u is just a single scalar value.
    bottom = p.k + p.N * p.beta * x(2) + p.b * x(2);
    
    % Compute mu = effect of therapy, mitigated by some resistance factors
    mu = top/bottom;
    
    % Compute population rate of change dx/dt
    dx(1) = x(1).*(p.r*((K-x(1))./K)-mu);

    % Compute strategy value rate of change du/dt
    dx(2) = p.s*(-(p.r*x(1).*x(2))/(p.sig*p.Kmax*exp(x(2).^2./(2*p.sig^2))) ...
        + (top*p.b)/(p.k+p.N*p.beta*x(2)+p.b*x(2))^2);
end

