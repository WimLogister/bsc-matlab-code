function h = get_fitness_handle( sys_input, time_points, rk_steps )
% Calling the outer function get_fitness_handle sets up everything to solve
% the model given system input (initial values and time range), the model
% parameters and number of integration steps. It returns a handle to anonymous
% function evaluate_fitness, which can subsequently be called by the optimization
% toolbox to evaluate the fitness for a given treatment schedule.

global counter
global soltab
global params

counter = 0;

h = @evaluate_fitness;
    function fit = evaluate_fitness( treat_sched )
    % Given a treatment schedule, this anonymous function solves the
    % ODE system and evaluates its fitness (average population over total time
    % range at the moment. This is the function that the optimization toolbox
    % uses to optimize the model.
        
        % Set up a treatment handle th by providing the number of time points
        % and a treatment schedule. This handle can then be used by RK to check
        % the treatment dosage at a given time t
        th = get_treatment_handle( time_points, treat_sched );
        
        t=0;
        dt=sys_input.tmax/rk_steps;
        
        sol=[sys_input.x0, sys_input.u0];
        soltab(1,:) = [t, sol, th(t)];
        
        % soltab consists of 4 columns:
        % 1. vector of time points t
        % 2. vector of population size at time t
        % 3. vector of resistance strategy at time t
        % 4. vector of intensity of chemotherapy at time t
        
        % Solve the system using RK4 integration for given number of steps
        for i=2:rk_steps
            sol=rk(t,sol,dt,th);
            t=t+dt;
            soltab(i,:) = [t,sol,th(t)];
        end
        
        % A counter that outputs the number of thousands of times the fitness
        % function is evaluated
        counter = counter + 1;
        if mod(counter,1000) == 0
            disp(counter);
        end
        
        % The fitness given the current model setup and treatment schedule,
        % represented by the average number of cancer cells over the total
        % simulation time
        fit = sum(soltab(:,2));
    end
end
