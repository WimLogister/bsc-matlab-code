%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. Declare system parameters and input %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    optimize = 1; % 0 For regular solving, > 0 for optimization
    create_plot = 0; % 0 For suppressing plot, > 0 for showing plot
    save_data = 1; % 0 For discarding plot, > 0 for saving plot to disk
    outer_loop = 1:numel(Ns); % Outer loop controlling N
    mid_loop = 1:numel(alphas_betas); % Middle loop controls how often we increment parameter of interest
    treatnum = 100;
    save_results = 0; % > 0 To save optimized treatment schedule, 0 to discard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = outer_loop % Outer loop controlling different N
    currN = Ns{i};
    params.N = currN.Nfun;
    
    for j = mid_loop % Middle loop controlling alpha and beta values
        eff = alphas_betas{j};
        params.alpha = eff.alpha;
        params.beta = eff.beta;
        
        for k = 1:1 % Inner loop controlling control points

            soltab = [];
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
                res=fmincon(treat,m0,A,b,[],[],zeros(1,numel(m0)),mmax*ones(1,numel(m0)),[],options); 
            else
                res=treat(m0cons);
                cons_tag = '_CONS';
            end

            muX=mean(soltab(:,2)); % Mean population density
            muU=mean(soltab(:,3)); % Mean resistance strategy
            
            if save_data > 0
                % Write optimized data to file
                filename=sprintf('%s_%s%s.m',eff.name,currN.name,cons_tag);
                dlmwrite(filename, soltab);
                
                % Write constant treatment data to file (ugly but it should
                % work)
                soltab = [];
                res=treat(m0cons);
                cons_tag = '_CONS';
                filename=sprintf('%s_%s%s.m',eff.name,currN.name,cons_tag);
                dlmwrite(filename, soltab);
            end
            
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

                sum_res=sum(res)
            end
        end
    end
end
