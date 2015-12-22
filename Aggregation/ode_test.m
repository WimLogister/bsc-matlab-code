clear variables

% System constants
rval=1.7; % Cancer growth rate
sigval=10; % Penalty to total pop. for increased resistance
Kmaxval=100; % Maximum carrying capacity
kval=0.005; % Cells' de novo resistance to therapy
bval=0.05; % Effectiveness of resistance
mval=0.1; % Chemotherapy dosage (paper says 0.1 for monotherapy)
sval=0.1; % Evolutionary speed

% Aggregation parameters
alphaval=1; % Power to determine the type of aggregation effect
betaval=0.6; % Scaling factor for other tumor cells' resistance
Nval=5; % Neighbourhood size

% Optimization parameters
mmax = 0.5;

% Model different aggregation effects by setting parameters as follows:
% 1. Dilution: alpha = beta = 0
% 2. Group detoxification: alpha = 1, beta > 0
% 3. Danger in numbers: alpha = 1.5, beta = 0
% 4. Group sellout: beta < 0

global params

params = struct('r',rval,'sig',sigval,'Kmax',Kmaxval,'k',kval,'b',bval,...
    'm',mval,'s',sval,'alpha',alphaval,'beta',betaval,'N',Nval);

dilution = struct('name','dilution','alpha',0,'beta',0,'rksteps',2000);
group_detox = struct('name','group detoxification','alpha',1,'beta',0.1,'rksteps',30000);
danger = struct('name','danger in numbers','alpha',1.4,'beta',0,'rksteps',20000);
sellout = struct('name','group sellout','alpha',1,'beta',-0.005,'rksteps',5000);
switchover = struct('name','switchover_LOWER_M','alpha',0,'beta',-0.005,'rksteps',5000);
best_case = struct('name','best_case','alpha',0,'beta',0.1,'rksteps',5000);
worst_case = struct('name','worst_case','alpha',1.4,'beta',-0.005,'rksteps',5000);

alphas_betas = {dilution group_detox danger sellout switchover best_case};
%alphas_betas = {dilution};
%alphas_betas = {group_detox};
%alphas_betas = {danger};
%alphas_betas = {sellout};
%alphas_betas = {switchover};
%alphas_betas = {best_case};
%alphas_betas = {worst_case};

N1 = struct('name','N=1','Nfun',@(x)1);
N5 = struct('name','N=5','Nfun',@(x)5);
N10 = struct('name','N=10','Nfun',@(x)10);
Nx = struct('name','N=1+x_10','Nfun',@(x)1+x/10);
Nxx = struct('name','N=x','Nfun',@(x)x);
Nx2 = struct('name','N=x_2','Nfun',@(x)x/5);
Ns = {Nx};

therapy;