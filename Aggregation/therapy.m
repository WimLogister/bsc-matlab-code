%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. Declare system parameters and input %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System constants
rval=0.1; % Cancer growth rate
sigval=1; % Penalty to total pop. for increased resistance
Kmaxval=100; % Maximum carrying capacity
kval=0.1; % Cells' de novo resistance to therapy
bval=5; % Effectiveness of resistance
mval=0.1; % Chemotherapy dosage (paper says 0.1 for monotherapy)
sval=0.001; % Evolutionary speed

% Aggregation parameters
alphaval=1; % Power to determine the type of aggregation effect
betaval=0.6; % Scaling factor for other tumor cells' resistance
Nval=5; % Neighbourhood size

% Model different aggregation effects by setting parameters as follows:
% 1. Dilution: alpha = beta = 0
% 2. Group detoxification: alpha = 1, beta > 0
% 3. Danger in numbers: alpha = 1.5, beta = 0
% 4. Group sellout: beta < 0

params = struct('r',rval,'sig',sigval,'Kmax',Kmaxval,'k',kval,'b',bval,...
    'm',mval,'s',sval,'alpha',alphaval,'beta',betaval,'N',Nval);

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
tmax=10000; % Total simulation time
steps=10000; % Number of integration steps used in ODE solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 2. Sample treatment strategies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Different treatment strategies are represented by anonymous functions of
% time that should be passed as an argument to the derivative function
% dosedyn. dosedyn Is the function that is passed to the ODE solver.

% Uniform constant treatment
treat0 = @(t) 1;

% Constant treatment with low intensity phase followed by high intensity phase
treat1 = @(t) 0.25*(t >= 0 & t <= 2500) + ...
    1.25*( t >= 2501 & t <= 10000);

% Constant treatment of increasing intensity
treat2 = @(t)0.4*(t >= 0 & t < 2500) + 0.8*(t >= 2500 & t < 5000) ...
    + 1.2*(t >= 5000 & t < 7500) + 1.6*(t >= 7501 & t < 10000);

% Linearly increasing treatment strategy
treat3 = @(t) t*2/tmax;

% Quadratically increasing treatment strategy
treat4 = @(t) t^2*(3/tmax^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 3. Solve system and display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the ODE system
[T,X] = ode45(@(t,x) dosedyn(t,x,treat0(t),params),0:tmax/steps:tmax,[tumorIni stratIni]);

muX=mean(X(:,1)) % Average tumor population size
muU=mean(X(:,2)) % Average resistance strategy value

x_label = sprintf('mu_{X} = %.3f',muX);
u_label = sprintf('mu_{u} = %.3f',muU);

% Plot results
figure(1)

% Population subplot
subplot(211), plot(T,X(:,1),'r'), line([0 tmax], [muX muX], 'Color', 'k')
legend('population density'),axis([0 tmax 0 100])
text(8000,muX+10,x_label)

% Resistance strategy subplot
subplot(212), plot(T,X(:,2),'b'), line([0 tmax], [muU muU], 'Color', 'k')
legend('resistance value', 'Location', 'northwest')
text(8000,muU+max(X(:,2))/10,u_label)