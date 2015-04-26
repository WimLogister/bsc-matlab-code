function h = get_fitness_handle( x0, u0, params )
h = @evaluate_fitness;
    function fit evaluate_fitness( treat )
        [T,X] = ode45(@(t,x) dosedyn(t,x,treat,params),linspace(0,10000,100),[x0 u0]);
        fit = sum(X(1,:));
    end
end
