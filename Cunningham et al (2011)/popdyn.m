function [ dx ] = popdyn( t,x )
    % Population and strategy dynamics for a single scalar phenotype strategy
    global r k be bp m sig Kmax schedule index treatment
    
    % This will hold the change rates of population and strategy dynamics
    dx=zeros(2,1);
    
    % Compute penalty to carrying capacity due to evolved resistance
    K=Kfun(x(2,:));
    
    if t > schedule(index)
        index=index+1;
        treatment=~treatment;
    end

    % Total sum of population
    X=sum(x(1,:));
    
    dosage=m*treatment;
    
    % Compute mu = effect of therapy, mitigated by some resistance factors
    mu=dosage./(k+be+bp*x(2,:));
    
    % Compute population change rate
    dx(1) = x(1,:).*(r*((K-X)./K)-mu);
    % Set evolutionary speed
    s = 0.2;
    % Compute strategy value change rate
    dx(2) = s*(-(r*X.*x(2,:))/(sig*Kmax*exp(x(2,:).^2./(2*sig^2))) +(dosage*bp)/(k+be+bp*x(2,:))^2);
end
