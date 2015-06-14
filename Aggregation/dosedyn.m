function [ dx ] = dosedyn( t,x,treat )
% Population and strategy dynamics for a single scalar phenotype strategy
% according to Brown et al (2015).
% x(1) is the population, x(2) is the resistance value, treat is the
% treatment intensity, p is a struct array storing system parameters

global params

    % Stores dx/dt and du/dt
    dx=zeros(1,2);
    
    % Series of checks to keep population size and resistance amount within
    % reasonable bounds
    if x(1) > 100
        x(1) = 100;
    end
    
    if x(1) < 0
        x(1) = eps;
    end
    
    if x(2) < 0
        x(2) = 0;
    end
    
    % Compute carrying capacity K
    K=params.Kmax.*exp((-x(2).^2)./(2*params.sig^2));
    
    % K==0 leads to numerical problems (NaN)
    if K == 0
        K = eps;
    end
    
    % Compute numerator of treatment efficacy mu
    top = (treat*params.N(x(1))^params.alpha)/params.N(x(1));
    
    % Compute denominator of mu
    % Note: since all models we consider only use a single scalar strategy, u = v
    % and u is just a single scalar value.
    bottom = params.k + params.N(x(1)) * params.beta * x(2) + params.b * x(2);
    
    % Compute mu = effect of therapy, mitigated by some resistance factors
    mu = top/bottom;
    
    % Compute population rate of change dx/dt
    dx(1) = x(1).*(params.r*((K-x(1))./K)-mu);

    % Compute strategy value rate of change du/dt
    dx(2) = params.s*(-(params.r*x(1).*x(2))/(params.sig*params.Kmax*exp(x(2).^2./(2*params.sig^2))) ...
        + (top*params.b)/(params.k+params.N(x(1))*params.beta*x(2)+params.b*x(2))^2);
    
    if isnan(mu) || isnan(x(1)) || isnan(dx(1)) || isnan(dx(2))
        exp(1);
    end
end

