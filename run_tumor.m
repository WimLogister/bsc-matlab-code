filename='basic_no_resistance.csv';
pars=csvread(filename);

% Declare system parameters as global
rval=0.1; % Cancer growth rate
sigval=5; % Penalty to total pop. for increased resistance
kval=0.1; % Baseline (de novo) resistance (not considering evolved resistance)
beval=0; % Resistance due to environmental factors
bpval=5; % Effectiveness of resistance strategy
mval=0.1/5; % Chemotherapy dosage (paper says 0.1 for monotherapy)
Kmaxval=100; % Maximum carrying capacity
uval=0; % Evolved resistance
setGlobalParams(rval,sigval,kval,beval,bpval,mval,Kmaxval,uval);

% Declare system input variables
tumorIni=100;

% Solve the ODE system
[T,X] = ode45(@popdyn,[0 250],tumorIni);

% Plot results
figure(1)
plot(T,X,'r')
title('Response of tumor population density to chemotherapy, no evolved resistance')