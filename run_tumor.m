% Declare system parameters as global
rval=0; % Cancer growth rate
sigval=0; % Penalty to total pop. for increased resistance
kval=0; % Baseline resistance (not considering evolved resistance)
beval=0; % Resistance due to environmental factors
bpval=0; % Effectiveness of resistance strategy
mval=0; % Chemotherapy dosage
Kmaxval=100; % Maximum carrying capacity
vval=0; % Evolved resistance
setGlobalParams(rval,sigval,kval,beval,bpval,mval,Kmaxval,vval);

% Declare system input variables
tumorIni=KMax;

% Solve the ODE system
[T,X] = ode45(@tumor,[0 250],tumorIni);

% Plot results
figure(1)
plot(T,X,'r')
title('Response of tumor population density to chemotherapy, no evolved resistance')