figure()
scatter(xM,yM,'r','^')
hold on
scatter(xK,yK)
%p = GreedyPilotSet;
% for i=1:20
%     matrix = find(pilot == p(:,i));
%     idx(i) = fix(matrix(1,1)/10)+1 ;
% end
idx = pilot_genetic;

for i=1:length(xK) 
    text(xK(i)+10,yK(i)+10,int2str(idx(i)))
end
X = [[xK;xM] [yK;yM]];
numcul = 4;
%rng('default') % For reproducibility
opts = statset('Display','final');
[id,C] = kmeans(X,numcul,'Distance','cityblock',...
    'Replicates',5,'Options',opts);
gscatter(X(:,1),X(:,2),id)
title(strcat(int2str(K),' Users',int2str(M),' APs Average Downlink Rate'))
xlabel('Distance (1000m)');
ylabel('Distance (1000m)');
idx1  = zeros(1,K);
for i=1:numcul
    l =[];
    for j=1:K
       if(id(j) == i)
          l = [l j];  
       end
       a = length(l);
       if(a<=10)
          matrix =  randperm(10,a);
       else
          m =  randperm(10,10);
          matrix = [m randperm(10,a-10)]
       end
       for k=1:a
           idx1(l(k)) = matrix(k);
       end
    end
end
Rate_Greedy_Avrg1 = [];
PilotSet_cluster = [];
pilotcontamination_greedy1 = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet_cluster = [PilotSet_cluster pilot(:,idx1(j))];
end
PilotSet = PilotSet_cluster;
mean(functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations))
for i=1:iteration
[PilotContamination_avgr, subMConta1] = fitness_contamination_forgreedy(PilotSet,Beta,K,M,nbrOfRealizations);
[C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
Rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
[GreedyPilotSet PilotSet] = functionGreedyPilotAssignment(M,K,Rate,Beta,pilot,PilotSet,nbrOfRealizations);
Rate_Greedy = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, GreedyPilotSet, nbrOfRealizations);
Rate_Greedy_Avrg1 = [Rate_Greedy_Avrg1 mean(Rate_Greedy)];
[PilotContamination_avgr1, subMConta1] = fitness_contamination_forgreedy(GreedyPilotSet,Beta,K,M,nbrOfRealizations);
pilotcontamination_greedy1 = [pilotcontamination_greedy1 PilotContamination_avgr];
end
for i=1:length(xK) 
    text(xK(i)+10,yK(i)-20,int2str(idx1(i)))
end



























