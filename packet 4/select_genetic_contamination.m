function [parent] = select_genetic_contamination(pop,pilot,Beta,K,M,nbrOfRealizations)
fitness =[];
l = size(pop);
for i=1:l(1)
    fit = fitness_contamination(pop(i,:),pilot,Beta,K,M,nbrOfRealizations);
    fitness = [fitness fit];
end
parent = [];
for i=1:l(1)/2
    index_min = find(fitness==min(fitness));
    parent = [parent;pop(index_min(1),:)];
    fitness(index_min) = 999999;
end

