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
Nval=5; % Neighbourhood size

% Dilution: alpha = beta = 0
% Group detoxification: alpha = 1, beta > 0
% Danger in numbers: alpha = 1.5, beta = 0
% Group sellout: beta < 0

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
tmax=10000;
steps=10000;

% Lengths of treatment and rest periods
treat_length=10;
rest_length=2000;

% Treatment bookkeeping variables
treatmentval=1;
indexval=1;
schedule=treat_sched(steps,treat_length,rest_length);
schedule=[250 1000 1250 2000 2250 3000 3250 4000 4250 5000 5250 6000 6250 7000 7250 8000 8250 9000 9250 10000];
%schedule=[500 2500 3000 5000 5500 7500 8000 9000 9500 10000];
%schedule=[1 7500 10000];
%schedule=[tmax];

% Declare constants and parameters as global
setGlobalParams(rval,sigval,alphaval,Nval,kval,bval,betaval,mval,Kmaxval,...
    sval,indexval,treatmentval,schedule);

% Solve the ODE system
[T,X] = ode45(@(t,x) aggdyn(t,x),0:tmax/steps:tmax,[tumorIni stratIni]);
%[T,X] = ode45(@(t,x) aggdyn(t,x),[0 tmax],[tumorIni stratIni]);
Xs=X(:,1);

mu=mean(X(:,1))

% Plot results
figure(1)
plot(T,X(:,1),'r')
axis([0 tmax 0 100]);
title('population density')
figure(2)
plot(T,X(:,2),'b');
title('average resistance strategy')