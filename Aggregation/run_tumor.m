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

% Declare constants and parameters as global
setGlobalParams(rval,sigval,alphaval,Nval,kval,bval,betaval,mval,Kmaxval,sval);

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
t=10000;

% Solve the ODE system
[T,X] = ode45(@aggdyn,[0 t],[tumorIni stratIni]);
Xs=X(:,1);

% Dilution: alpha = beta = 0
% Group detoxification: alpha = 1, beta > 0
% Danger in numbers: alpha = 1.5, beta = 0
% Group sellout: beta < 0


% Plot results
figure(1)
plot(T,Xs,'r')
axis([0 t 0 100]);
title('Response of tumor population density to chemotherapy, with evolved resistance')