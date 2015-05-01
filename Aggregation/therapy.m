%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. Declare system parameters and input %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

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

global params

params = struct('r',rval,'sig',sigval,'Kmax',Kmaxval,'k',kval,'b',bval,...
    'm',mval,'s',sval,'alpha',alphaval,'beta',betaval,'N',Nval);

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

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
tmax=10000; % Total simulation time

treatnum=20; % Number of treatment / control points
m0=zeros(1,treatnum); % Initial treatment guess for optimizer
%m0=0.1+m0;
global counter
global soltab

% Save system input in structure
system_input=struct('x0',tumorIni,'u0',stratIni,'tmax',tmax);
T=linspace(0,9900,treatnum); % Vector holding time points
rk_timesteps=5000; % Number of Runge-Kutte ODE integration steps

h = get_fitness_handle(system_input,T,rk_timesteps);

options=optimset('Maxiter',1000,'DiffMinChange',1e-6);
res=fmincon(h,m0,[],[],[],[],zeros(1,numel(m0)),0.1*ones(1,numel(m0)),[],options);

figure(1)

plot(soltab(:,1),soltab(:,2))
