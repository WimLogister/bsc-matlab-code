% Declare system parameters as global
rval=0.1; % Cancer growth rate
sigval=5; % Penalty to total pop. for increased resistance
kval=0.1; % Baseline (de novo) resistance (not considering evolved resistance)
beval=0; % Resistance due to environmental factors
bpval=5; % Effectiveness of resistance strategy
mval=0.1; % Chemotherapy dosage (paper says 0.1 for monotherapy)
Kmaxval=100; % Maximum carrying capacity

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
tmax=250;

% 100 time steps of treatment, 100 time steps of rest
treat_length=50;
rest_length=50;

% Treatment bookkeeping variables
treatmentval=1;
indexval=1;
schedule=treat_sched(tmax,treat_length,rest_length);
%schedule=[tmax];

setGlobalParams(rval,sigval,kval,beval,bpval,mval,Kmaxval,indexval,treatmentval,schedule);

% Solve the ODE system
[T,X] = ode45(@popdyn,[0 tmax],[tumorIni stratIni]);
Xs=X(:,2);

% Plot results
figure(1)
plot(T,Xs,'r')
%axis([0 250 0 100]);
axis tight
title('Response of tumor population density to chemotherapy, with evolved resistance')