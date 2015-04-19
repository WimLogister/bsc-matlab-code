function [ schedule ] = treat_sched( length, treat_length, rest_length )
    index=0;
    treat=1; % 1 for treatment period, 0 for rest period
    schedule=[];
    while index < length
        if treat > 0
            index=index+treat_length;
        else
            index=index+rest_length;
        end
        schedule=[schedule index];
        treat=~treat;
    end
end
