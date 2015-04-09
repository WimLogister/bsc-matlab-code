% Declare system parameters as global
rval=0.1; % Cancer growth rate
sigval=1; % Penalty to total pop. for increased resistance
Kmaxval=100; % Maximum carrying capacity
alphaval=0; % Power to determine the type of aggregation effect
Nval=10; % Neighbourhood size
kval=0.1; % Cells' de novo resistance to therapy
bval=5; % Effectiveness of resistance
betaval=0; % Scaling factor for other tumor cells' resistance
mval=0.1; % Chemotherapy dosage (paper says 0.1 for monotherapy)


setGlobalParams(rval,sigval,alphaval,Nval,kval,bval,betaval,mval,Kmaxval);

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value

% Solve the ODE system
[T,X] = ode45(@aggdyn,[0 1000],[tumorIni stratIni]);
Xs=X(:,1);

% Plot results
figure(1)
plot(T,Xs,'r')
axis([0 1000 0 100]);
title('Response of tumor population density to chemotherapy, with evolved resistance')