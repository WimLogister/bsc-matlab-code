function h = get_fitness_handle( sys_input, time_points,  rk_steps )
% Calling the outer function get_fitness_handle sets up everything to solve
% the model given system input (initial values and time range), the model
% parameters and number of integration steps. It returns a handle to anonymous
% function evaluate_fitness, which can subsequently be called by the optimization
% toolbox to evaluate the fitness for a given treatment schedule.

% What does this function need? rk needs to be called somewhere

global counter
global soltab

counter = 0;

h = @evaluate_fitness;
    function fit = evaluate_fitness( treat_sched )
    % Given a treatment schedule, this anonymous function solves the
    % ODE system and evaluates its fitness. This is the function that the
    % optimization toolbox will use to optimize the model.
        
        % Treatment handle th can be used by RK to check the treatment
        % dosage at a given time t
        th = get_treatment_handle( time_points, treat_sched );
        
        t=0;
        dt=sys_input.tmax/rk_steps;
        
        sol=[sys_input.x0, sys_input.u0];
        soltab(1,:) = [t, sol, th(t)];
        
        for i=2:rk_steps
            sol=rk(t,sol,dt,th);
            t=t+dt;
            soltab(i,:) = [t,sol,th(t)];
        end
        
        counter = counter + 1;
        if mod(counter,100) == 0
            disp(counter);
        end
        
        fit = sum(soltab(:,2));
    end
end
