%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 2. Solve system and display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numel(Ns) % Outer loop that loops over different values for N
    currN = Ns{i}; % Store current N structure for easy identification
    params.N = currN.Nfun; % In parameter structure, set currN as the N to use
    
    for j = 1:numel(effects) % Middle loop that loops over different effects
        eff = effects{j};
        params.alpha = eff.alpha; % Store value of alpha in parameter struct
        params.beta = eff.beta; % Store value of beta in parameter struct

        soltab = [];
        m0=zeros(1,treatnum); % Initial treatment guess for optimizer

        T=linspace(0,tmax-tmax/treatnum,treatnum); % Vector holding time points

        % 'fitness' is a handle to the fitness evaluation function that
        % fmincon will call during optimization.
        % We pass the system input, the number of time points and
        % number of Runge-Kutta integration steps, which will be stored
        % in the 'fitness' handle.
        fitness = get_fitness_handle(system_input,T,eff.rksteps);

        m0cons=[params.m*ones(1,treatnum)]; % Constant treatment schedule
        %m0cons=zeros(1,treatnum);
        %m0cons(5:25)=0.5;

        % A and b are a system of equations that is passed to fmincon to
        % make sure that the total treatment sums up to the appropriate amount
        A=ones(1,treatnum);
        b=params.m*treatnum;

        if optimize > 0
            % Optimize treatment
            options=optimset('Maxiter',1000,'DiffMinChange',1e-12);
            res=fmincon(fitness,m0,A,b,[],[],zeros(1,numel(m0)),mmax*ones(1,numel(m0)),[],options); 
        else
            % Apply constant treatment
            res=fitness(m0cons);
        end

        muX=mean(soltab(:,2)); % Mean population density
        muU=mean(soltab(:,3)); % Mean resistance strategy

        if save_data > 0
            if optimize > 0
                % Write optimized data to file
                filename=sprintf('%s_%s_%s.m',eff.name,currN.name,ID_tag);
                dlmwrite(filename, soltab);
            end

            % Write constant treatment data to file (ugly but it works)
            soltab = [];
            res=fitness(m0cons);
            cons_tag = '_CONS';
            filename=sprintf('%s_%s_CONS_%s.m',eff.name,currN.name,ID_tag);
            dlmwrite(filename, soltab);
        end

        if make_plot > 0
            create_plot; % Plot treatment schedule
        end
    end
end
