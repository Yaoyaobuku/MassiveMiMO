function [parent] = select_genetic(pop,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations)
fitness =[];
l = size(pop);
for i=1:l(1)
    fit = fitness_rate(pop(i,:),pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    fitness = [fitness fit];
end
parent = [];
for i=1:l(1)/2
    index_max = find(fitness==max(fitness));
    parent = [parent;pop(index_max(1),:)];
    fitness(index_max) = -999999;
end
