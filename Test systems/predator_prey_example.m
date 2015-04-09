% Declare system parameters a, b, c and d as global variables
setGlobalParams(0.4,0.001,0.001,0.9);

% Declare input variables
Ntot=1000;
p=0.4;
RabbitIni=(1-p)*Ntot;
FoxIni=p*Ntot;

% Solve the ODE system
[T,Y]=ode45(@rabbitfox,[0 200],[RabbitIni FoxIni]);

% Plot results
figure(1)
plot(T,Y(:,1),'r')
hold on
plot(T,Y(:,2),'b')
legend('Rabbits','Foxes')
title('Predator-Prey model ODE')