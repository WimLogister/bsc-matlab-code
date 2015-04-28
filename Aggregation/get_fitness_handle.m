function h = get_fitness_handle( sys_input, model_params, time_points,  rk_steps )
% Calling the outer function get_fitness_handle sets up everything to solve
% the model given system input (initial values and time range), the model
% parameters and number of integration steps. It returns a handle to anonymous
% function evaluate_fitness, which can subsequently be called by the optimization
% toolbox to evaluate the fitness for a given treatment schedule.

% What does this function need? rk needs to be called somewhere

h = @evaluate_fitness;
    function fit = evaluate_fitness( treat_sched )
    % Given a treatment schedule, this anonymous function solves the
    % ODE system and evaluates its fitness. This is the function that the
    % optimization toolbox will use to optimize the model.
        
        % Treatment handle th can be used by RK to check the treatment
        % dosage at a given time t
        th = get_treatment_handle( time_points, treat_sched );
        
        % Then runge-kutta needs
        % to be passed the treatment handle, so that it can make the
        % necessary calls to treat(t)
        
        t=0;
        dt=sys_input.tmax/rk_steps;
        
        for i=1:rk_steps
            sol=rk(t,sol,dt,th);
            t=t+dt;
            % need to implement some kind of table to store the solutions
            soltab;
        end
        
        %[P,X] = ode45(@(t,x) dosedyn(ttable(t),x,treat,params),L,[x0 u0]);
        %fit = X;
        %fit = sum(X(1,:));
    end
end
