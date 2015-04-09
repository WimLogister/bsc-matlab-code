function [ G ] = Gfun( v )
    global r alpha N k b beta m
   
    x = 100;
    K=Kfun(v);
    
    % Numerator of mu
    top = m*N^alpha/N;
    
    % Denominator of mu
    % Note: since all models we consider only use a single scalar strategy, u = v
    % and u is just a single scalar value.
    bottom = k + N * beta * v + b * v;
    
    % Compute mu = effect of therapy, mitigated by some resistance factors
    mu = top/bottom;
    
    G = r*(K-x)/K - mu;
end

