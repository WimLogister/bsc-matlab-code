% Declare system parameters as global
rval=0.1; % Cancer growth rate
sigval=5; % Penalty to total pop. for increased resistance
kval=0.1; % Baseline (de novo) resistance (not considering evolved resistance)
beval=1; % Resistance due to environmental factors
bpval=5; % Effectiveness of resistance strategy
mval=0.08; % Chemotherapy dosage (paper says 0.1 for monotherapy)
Kmaxval=100; % Maximum carrying capacity
vval=0; % Evolved resistance
setGlobalParams(rval,sigval,kval,beval,bpval,mval,Kmaxval,vval);

% Declare system input variables
tumorIni=100;

% Solve the ODE system
[T,X] = ode45(@tumor,[0 250],tumorIni);

% Plot results
figure(1)
plot(T,X,'r')
title('Response of tumor population density to chemotherapy, no evolved resistance')