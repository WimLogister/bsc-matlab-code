function h = get_treatment_handle( t_i, m_i )
% Based on input arguments, get_treatment_handle sets up a treatment schedule
% and returns a handle h for the anonymous inner function treat. treat can
% then be used to poll the treatment schedule by simply passing a time t.

% T is a vector of time points t_i, m is a vector of intensities of the treatment
% in the intervals between the t_i in T.
h = @treat;
    function dose = treat( t )
        dose = 0;
        if t <= t_i(1)
            dose = m_i(1);
        else if t >= t_i(end)
            dose = m_i(end);
        else
            i = 1;
            while t > t_i(i)
                i=i+1;
            end
            dose = m_i(i-1);
            end
        end
    end

end

