%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. Declare system parameters and input %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

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

% Model different aggregation effects by setting parameters as follows:
% 1. Dilution: alpha = beta = 0
% 2. Group detoxification: alpha = 1, beta > 0
% 3. Danger in numbers: alpha = 1.5, beta = 0
% 4. Group sellout: beta < 0

global params

params = struct('r',rval,'sig',sigval,'Kmax',Kmaxval,'k',kval,'b',bval,...
    'm',mval,'s',sval,'alpha',alphaval,'beta',betaval,'N',Nval);

dilution = struct('name','dilution','alpha',0,'beta',0,'rksteps',10000);
group_detox = struct('name','group detoxification','alpha',1,'beta',0.1,'rksteps',50000);
danger = struct('name','danger in numbers','alpha',1.3,'beta',0,'rksteps',10000);
sellout = struct('name','group sellout','alpha',1,'beta',-0.005,'rksteps',50000);
switchover = struct('name','switchover','alpha',0,'beta',-0.005,'rksteps',50000);
best_case = struct('name','best_case','alpha',0,'beta',0.1,'rksteps',50000);
worst_case = struct('name','worst_case','alpha',1.3,'beta',-0.005,'rksteps',50000);

%alphas_betas = {dilution};

%alphas_betas = {dilution group_detox danger sellout switchover best_case worst_case};
alphas_betas = {dilution danger};

N1 = struct('name','N=1','Nfun',@(x)1);
N5 = struct('name','N=5','Nfun',@(x)5);
N10 = struct('name','N=10','Nfun',@(x)10);
Nx = struct('name','N=1+x_10','Nfun',@(x)1+x/10);
Nxx = struct('name','N=x','Nfun',@(x)x);
Ns = {Nxx};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 2. Solve system and display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare system input variables
tumorIni=100; % Initial cancer cell population size
stratIni=0.0; % Initial phenotypic strategy (resistance) value
tmax=100; % Total simulation time

global soltab % Used to store solutions to differential equations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation control variables                                  %
    optimize = 0; % 0 For regular solving, > 0 for optimization
    show_plot = 1; % 0 For suppressing plot, > 0 for showing plot
    save_plot = 0; % 0 For discarding plot, > 0 for saving plot to disk
    outer_loop = 1:numel(Ns); % Outer loop controlling N
    mid_loop = 1:numel(alphas_betas); % Middle loop controls how often we increment parameter of interest
    inner_loop = 0.5; % Inner loop controls how often we increment # of control pts
    treatnum = 100;
    save_results = 0; % > 0 To save optimized treatment schedule, 0 to discard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

create_plot = 0;
if show_plot > 0 || save_plot > 0
    create_plot = 1;
end

for i = outer_loop % Outer loop controlling different N
    currN = Ns{i};
    params.N = currN.Nfun;
    
    for j = mid_loop % Middle loop controlling alpha and beta values
        eff = alphas_betas{j};
        params.alpha = eff.alpha;
        params.beta = eff.beta;
        
        for mmax = inner_loop % Inner loop controlling control points

            m0=zeros(1,treatnum); % Initial treatment guess for optimizer

            % Save system input in structure
            system_input=struct('x0',tumorIni,'u0',stratIni,'tmax',tmax);
            T=linspace(0,tmax-tmax/treatnum,treatnum); % Vector holding time points
            rk_timesteps=eff.rksteps; % Number of Runge-Kutte ODE integration steps

            % treat is a function handle to the fitness function
            treat = get_fitness_handle(system_input,T,rk_timesteps);
            treat_cutoff = 100;

            % A and b are a system of equations that make sure that all the decision
            % variables sum to the appropriate amount
            A=ones(1,treatnum);
            b=params.m*treatnum;
            m0cons=[params.m*ones(1,treat_cutoff)];
            
            %m0cons=[params.m*ones(1,treat_cutoff) zeros(1,treatnum-treat_cutoff)];
            %A=[ones(1,treat_cutoff) zeros(1,treatnum-treat_cutoff);...
            %    zeros(1,treat_cutoff) ones(1,treatnum-treat_cutoff)];
            %b=[params.m*treat_cutoff; 0];

            cons_tag = '';
            if optimize > 0
                options=optimset('Maxiter',1000,'DiffMinChange',1e-6);
                %options=optimset('Maxiter',1000,'DiffMinChange',1e-12,'Display','iter-detailed');
                res=fmincon(treat,m0,A,b,[],[],zeros(1,numel(m0)),mmax*ones(1,numel(m0)),[],options); 
                if save_results > 0
                    res_filename = sprintf('%utmax%ucpts_%s',tmax,treatnum,eff.name);
                    save(res_filename,'res');
                end
            else
                res=treat(m0cons);
                cons_tag = '_CONS';
            end

            muX=mean(soltab(:,2)); % Mean population density
            muU=mean(soltab(:,3)); % Mean resistance strategy
            
            filename=sprintf('%s_%s%s.m',eff.name,currN.name,cons_tag);
            dlmwrite(filename, soltab);

            if create_plot > 0
                % Create plot showing results
                my_fig=figure;

                x_label = sprintf('mu_{X} = %.3f',muX); % Label for population mean
                u_label = sprintf('mu_{u} = %.3f',muU); % Label for resistance mean
                m_title = sprintf('Opt. treat sched., %u treatment periods, mmax=%.3f',treat_cutoff,mmax); % Title for treatment plot
                x_title = sprintf('Pop dty vs time, %s, alpha=%.1f, beta=%.1f, effect=%s',currN.name,params.alpha,params.beta,eff.name);

                % Population subplot
                subplot(311),plot(soltab(:,1),soltab(:,2),'r'), line([0 tmax], [muX muX], 'Color', 'k')
                title(x_title),axis([0 tmax 0 params.Kmax])
                text(tmax-tmax*0.2,muX+1,x_label)

                % Resistance strategy subplot
                subplot(313), plot(soltab(:,1),soltab(:,3),'b'), line([0 tmax], [muU muU], 'Color', 'k')
                title('Evolved resistance vs time'),axis tight
                text(tmax-tmax*0.2,muU+max(soltab(:,3))/10,u_label)

                % Optimized treatment regime subplot
                subplot(312), plot(soltab(:,1),soltab(:,4),'g'),
                title(m_title),axis([0 tmax 0 mmax*1.25])

                % Save plot to file
                if save_plot > 0
                    %baseFileName = sprintf('figure_%d.jpg',k);
                    baseFileName = sprintf('%s_%s.png',currN.name,eff.name);
                    fullFileName = fullfile('C:\Users\Wim\Documents\KE\Bsc Thesis\Code\Aggregation\26-5-2015', baseFileName);  
                    my_fig; % Activate the figure again.
                    exportfig(fullFileName); % Using export_fig instead of saveas.
                end
                sum_res=sum(res)
            end
        end
    end
end
