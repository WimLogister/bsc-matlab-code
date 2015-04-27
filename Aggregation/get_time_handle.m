function h = get_time_handle( table )
h = @get_time_at;
    function time = get_time_at( t )
        time = table(t);
    end
end
