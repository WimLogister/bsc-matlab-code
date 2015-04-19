% System constants
rval=0.1; % Cancer growth rate
sigval=1; % Penalty to total pop. for increased resistance
Kmaxval=100; % Maximum carrying capacity
kval=0.1; % Cells' de novo resistance to therapy
bval=5; % Effectiveness of resistance
mval=0.1; % Chemotherapy dosage (paper says 0.1 for monotherapy)
sval=0.001; % Evolutionary speed

% System parameters
alphaval=0; % Power to determine the type of aggregation effect
betaval=0; % Scaling factor for other tumor cells' resistance
Nval=1; % Neighbourhood size

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
tmax=10000;
steps=10000;

% 100 time steps of treatment, 100 time steps of rest
treat_length=2000;
rest_length=1000;

% Treatment bookkeeping variables
treatmentval=1;
indexval=1;
schedule=treat_sched(steps,treat_length,rest_length);
%schedule=[tmax];

% Declare constants and parameters as global
setGlobalParams(rval,sigval,alphaval,Nval,kval,bval,betaval,mval,Kmaxval,...
    sval,indexval,treatmentval,schedule);

% Solve the ODE system
%[T,X] = ode45(@(t,x) aggdyn(t,x),0:tmax/steps:tmax,[tumorIni stratIni]);
[T,X] = ode45(@(t,x) aggdyn(t,x),[0 tmax],[tumorIni stratIni]);
Xs=X(:,1);

mu=mean(X(:,1));

% Dilution: alpha = beta = 0
% Group detoxification: alpha = 1, beta > 0
% Danger in numbers: alpha = 1.5, beta = 0
% Group sellout: beta < 0


% Plot results
figure(1)
plot(T,Xs,'r')
axis([0 tmax 0 100]);
title('Response of tumor population density to chemotherapy, with evolved resistance')