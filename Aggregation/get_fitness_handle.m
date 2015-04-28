function h = get_fitness_handle( x0, u0, model_params, rk_params )
h = @evaluate_fitness;
    function fit = evaluate_fitness( treat )
        % Use RK to integrate ODE
        % RK needs: number of time steps
        % This function needs to create a time vector and treatment
        % vector to pass to get_treatment_handle. Then runge-kutta needs
        % to be passed the treatment handle, so that it can make the
        % necessary calls to treat(t)
        
        %[P,X] = ode45(@(t,x) dosedyn(ttable(t),x,treat,params),L,[x0 u0]);
        %fit = X;
        %fit = sum(X(1,:));
    end
end
