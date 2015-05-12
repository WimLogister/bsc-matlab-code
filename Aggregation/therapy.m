%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. Declare system parameters and input %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

load 'res.mat'

% System constants
rval=0.1; % Cancer growth rate
sigval=1; % Penalty to total pop. for increased resistance
Kmaxval=100; % Maximum carrying capacity
kval=0.1; % Cells' de novo resistance to therapy
bval=5; % Effectiveness of resistance
mval=0.1; % Chemotherapy dosage (paper says 0.1 for monotherapy)
sval=0.01; % Evolutionary speed

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
tmax=365; % Total simulation time

treatnum = tmax;
cntr = 1;

global soltab % Used to store solutions to differential equations

while treatnum < treatnum + 1
    
    m0=zeros(1,treatnum); % Initial treatment guess for optimizer
    mmax=0.2; % Maximum treatment amount
    normtreat=0.04; % "Normal" treatment amount, so the amount given when
    % giving constant treatment

    % Save system input in structure
    system_input=struct('x0',tumorIni,'u0',stratIni,'tmax',tmax);
    T=linspace(0,tmax-tmax/treatnum,treatnum); % Vector holding time points
    rk_timesteps=750; % Number of Runge-Kutte ODE integration steps

    % treat is a function handle to the fitness function
    treat = get_fitness_handle(system_input,T,rk_timesteps);

    % A and b are a system of equations that make sure that all the decision
    % variables sum to the appropriate amount
    A=ones(1,treatnum);
    b=normtreat*treatnum;

    optimize = 1; % 0 For regular solving, > 0 for optimization
    show_plot = 1; % 0 For not plotting, > 0 for plotting

    if optimize > 0
        options=optimset('Maxiter',1000,'DiffMinChange',1e-12,'Display','iter-detailed');
        res=fmincon(treat,m0,A,b,[],[],zeros(1,numel(m0)),mmax*ones(1,numel(m0)),[],options); 
    else
        res=treat(normtreat+m0);
    end

    muX=mean(soltab(:,2)); % Mean population density
    muU=mean(soltab(:,3)); % Mean resistance strategy

    output = struct('ctrlpts',treatnum,'muX',muX,'muU',muU,'optimal',res);

    out_filename = sprintf('optimize%u.dat',treatnum);

    struct2csv(output, out_filename);

    if show_plot > 0
        % Plot results
        figure(cntr)

        x_label = sprintf('mu_{X} = %.3f',muX); % Label for population mean
        u_label = sprintf('mu_{u} = %.3f',muU); % Label for resistance mean
        m_title = sprintf('Opt. treat sched., %u treatment periods',treatnum); % Title for treatment plot

        
        % Population subplot
        subplot(311),plot(soltab(:,1),soltab(:,2),'r'), line([0 tmax], [muX muX], 'Color', 'k')
        title('Population density vs time'),axis([0 tmax 0 100])
        text(tmax-tmax*0.2,muX+10,x_label)

        % Resistance strategy subplot
        subplot(313), plot(soltab(:,1),soltab(:,3),'b'), line([0 tmax], [muU muU], 'Color', 'k')
        title('Evolved resistance vs time'),axis tight
        text(tmax-tmax*0.2,muU+max(soltab(:,3))/10,u_label)
        
        % Optimized treatment regime subplot
        subplot(312), plot(soltab(:,1),soltab(:,4),'g'),
        title(m_title),axis tight
    end

    treatnum=treatnum*2; % Number of treatment / control points
    cntr=cntr+1;
end

% a = struct('muX',73.5,'muU',3.5,'numpoints',24); % Create structure
% struct2csv(a,'test.dat') % Write data
% A=importdata('test.dat') % Read data