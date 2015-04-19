function [ schedule ] = create_treat_sched( length, treat_length, rest_length )
    schedule=zeros(length,1);
    index=1;
    treatment=1; % 1 if treatment period, 0 if rest period, start with treatment
    while index <= length
        period_length=0;
        if treatment > 0
            period_length=treat_length;
        else
            period_length=rest_length;
        end
        if index + period_length > numel(schedule)
            schedule(index:end)=treatment;
        else
            schedule(index:index+period_length)=treatment;
        end
        treatment=~treatment;
        index=index+period_length;
    end
end
