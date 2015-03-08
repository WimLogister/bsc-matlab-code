function [ output_args ] = Gfun( u, x )
global r
% The G-function computes the fitness of a cell with a certain strategy
% based on the overall cell density and on the effect of the chemotherapy
% on cells with that particular strategy.
% Dus als input krijg je een sized-n vector x waarin de gehele populatie
% verdeeld is. Elke entry i geeft aan hoeveel cellen strategie i
% aanhangen. De G-functie returnt dan een vector g met de fitness voor
% elk type strategie. Vermenigvuldig x met g om de rate of change dx te
% berekenen voor elk type cell. Dit is dus wat je returnt in je ODE voor
% cell growth.

end

