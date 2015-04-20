function [ dx ] = aggdyn( t,x )
    % Population and strategy dynamics for a single scalar phenotype strategy
    global r sig alpha N k b beta m Kmax s schedule index treatment
    
    % This will hold the change rates of population and strategy dynamics
    dx=zeros(2,1);
    
    % Compute penalty to carrying capacity due to evolved resistance
    K=Kfun(x(2));
    
    if x(1) > 100
        x(1) = 100;
    end
    
    if x(1) < 0
        x(1) = 0;
    end
    
    if x(2) < 0
        x(2) = 0;
    end
    
    if K == 0
        K = eps;
    end
    
    if t > schedule(index)
        index=index+1;
        disp('index:');
        index
        treatment=~treatment;
    end
    
    dosage=treatment*m;
    
    % Numerator of mu
    top = (dosage*N^alpha)/N;
    
    % Denominator of mu
    % Note: since all models we consider only use a single scalar strategy, u = v
    % and u is just a single scalar value.
    bottom = k + N * beta * x(2) + b * x(2);
    
    % Compute mu = effect of therapy, mitigated by some resistance factors
    mu = top/bottom;
    
    % Compute population change rate
    dx(1) = x(1).*(r*((K-x(1))./K)-mu);
    if isnan(dx(1)) | isinf(dx(1))
        disp('NaN detected')
    end
    % Compute strategy value change rate
    dx(2) = s*(-(r*x(1).*x(2))/(sig*Kmax*exp(x(2).^2./(2*sig^2))) + (top*b)/(k+N*beta*x(2)+b*x(2))^2);
end
