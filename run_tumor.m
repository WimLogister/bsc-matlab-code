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
uval=5; % Evolved resistance
xval=100;
setGlobalParams(rval,sigval,kval,beval,bpval,mval,Kmaxval,uval,xval);

% Declare system input variables
tumorIni=100;

% Solve the ODE system
[T,U] = ode45(@stratdyn,[0 250],uval);

% Plot results
figure(1)
plot(T,U,'r')
title('Response of tumor population density to chemotherapy, no evolved resistance')