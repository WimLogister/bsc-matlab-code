function h = get_fitness_handle( x0, u0, params, T )
h = @evaluate_fitness;
ttable = get_time_handle(T);
L = linspace(1,100,100);
    function fit = evaluate_fitness( treat )
        [P,X] = ode45(@(t,x) dosedyn(ttable(t),x,treat,params),L,[x0 u0]);
        fit = X;
        %fit = sum(X(1,:));
    end
end
