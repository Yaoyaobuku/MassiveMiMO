clear;
tau_cf= 10;
K = 20 ; M=40; nbrOfRealizations = 1; D_sqr = 1000; population = 9;
taud_sc = 20; tauu_sc = 20; BW = 20e6; NF_dB = 9;
AVErhod_cf = 200; AVErhou_cf = 100; AVErhop_cf = 100; 
iteration = 30; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DistanceControl = 'Uni'; % Control two solutions creat uniformly distributed
% 'Halton' use Halton sequenc, and 'Uni' use makdedist Uniformly-distribution
ShadowingControl = 'uncorrelated'; % Control two shadowing correlation model: 'uncorrelated' or 'correlated'
PowerControl = 'No'; % Two Power Control Mode: 'No' = without Power Control / 'Yes' = Max-Min Power Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d_MK xM yM xK yK] = functionDistance(M, K, D_sqr, DistanceControl, nbrOfRealizations);
figure()
scatter(xM,yM,'^')
hold on
for i=1:length(xK) 
    scatter(xK(i),yK(i))
    text(xK(i)+10,yK(i)+10,int2str(i))
end
[Beta PL z_MK] = functionLargeScaleFading(d_MK, M, K, ShadowingControl, nbrOfRealizations);
 % Beta = ones(M, K, nbrOfRealizations); % beta_mk = 1
[NoisePower rhod_cf rhou_cf rhop_cf rhod_sc rhou_sc rhoup_sc rhodp_sc] = functionNormalizedTransmitSNRs(M, K, BW, NF_dB, AVErhod_cf, AVErhou_cf, AVErhop_cf);
[Hchannel Gchannel Wnoise] = functionGchannelGenerating(M, K, tau_cf, Beta, nbrOfRealizations);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pilot = [[1 0 0 0];[0 1 0 0 ];[0 0 1 0];[0 0 0 0]] %generate pilot
pilot = functionRandomPilotAssignment(tau_cf, tau_cf, nbrOfRealizations);
random = [randi([1 tau_cf],1,K)]; %random greedy
%random = functioncreatpop(K,tau_cf,1,Beta);
PilotSet = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet = [PilotSet pilot(:,random(j))];
end
[PilotContamination_avgr, subMConta1] = fitness_contamination(random,pilot,Beta,K,M,nbrOfRealizations);
[C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
Rate_start = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
Pilotcontamination_start =  pilot_contamination(PilotSet,Beta,K,M,nbrOfRealizations);
Rate_Greedy_Avrg = mean(Rate_start);
pilotcontamination_greedy = Pilotcontamination_start;
for i=1:iteration
[PilotContamination_avgr, subMConta1] = fitness_contamination_forgreedy(PilotSet,Beta,K,M,nbrOfRealizations);
[C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
Rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
[GreedyPilotSet PilotSet] = functionGreedyPilotAssignment(M,K,Rate,Beta,pilot,PilotSet,nbrOfRealizations);
Rate_Greedy = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, GreedyPilotSet, nbrOfRealizations);
Rate_Greedy_Avrg = [Rate_Greedy_Avrg mean(Rate_Greedy)];
[PilotContamination_avgr, subMConta1] = fitness_contamination_forgreedy(GreedyPilotSet,Beta,K,M,nbrOfRealizations);
pilotcontamination_greedy = [pilotcontamination_greedy PilotContamination_avgr];
end

% %genentic pilot assigment with pilot contamination
% pop = [];
% PilotSet = [];
% avr_pilotcontam_genetic = Pilotcontamination_start;
% Rate_genetic_pc = mean(Rate_start);
% pop = functioncreatpop(K,tau_cf,population,Beta);
% pop = [random;pop];
% for i=1:iteration
%     parent = select_genetic_contamination(pop,pilot,Beta,K,M,nbrOfRealizations);
%     child = cross_over_genetic(parent,tau_cf);
%     pop = [parent;child];
% fitness = [];
% l = size(pop);
% for i=1:l(1)
%     [fit a] = fitness_contamination(pop(i,:),pilot,Beta,K,M,nbrOfRealizations);
%     fitness = [fitness fit];
% end
% index_genetic = find(fitness == min(fitness));
% pilot_genetic_1 = pop(index_genetic(1),:);
% PilotSet1 = [];
% for j=1:K
%        %PilotSet caculates rate
%        PilotSet1 = [PilotSet1 pilot(:,pilot_genetic_1(j))];
% end
% [pc subMConta1] = fitness_contamination(pilot_genetic_1,pilot,Beta,K,M,nbrOfRealizations);
% avr_pilotcontam_genetic = [avr_pilotcontam_genetic pc];
% [C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet1, nbrOfRealizations);
% Rate_genetic = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet1, nbrOfRealizations);
% Rate_genetic_pc = [Rate_genetic_pc mean(Rate_genetic)];
% end
% %

%genentic pilot assigment with max average rate
pop = [];
PilotSet = [];
Rate_genetic_rate =mean(Rate_start);
rate_pilotcontam_genetic = Pilotcontamination_start;
pop = functioncreatpop(K,tau_cf,population,Beta);
pop = [random;pop];
parent = [];
for i=1:iteration
    parent = select_genetic(pop,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    child = cross_over_genetic(parent,tau_cf);
    pop = [parent;child];
fitness = [];
l = size(pop);
for i=1:l(1)
    fit = fitness_rate(pop(i,:),pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    fitness = [fitness fit];
end
index_genetic = find(fitness == max(fitness));
pilot_genetic = pop(index_genetic(1),:);
PilotSet1 = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet1 = [PilotSet1 pilot(:,pilot_genetic(j))];
end

[pc subMConta1] = fitness_contamination_forgreedy(PilotSet1,Beta,K,M,nbrOfRealizations)
[C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet1, nbrOfRealizations);
Rate_genetic_rateavrg = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet1, nbrOfRealizations);
Rate_genetic_rate = [Rate_genetic_rate mean(Rate_genetic_rateavrg)];
rate_pilotcontam_genetic = [rate_pilotcontam_genetic pc];
end

%K-Mean Cluster Algorithm
Rate_Greedy_Avrg1 = mean(Rate_start);
PilotSet_cluster = [];
pilotcontamination_greedy1 = Pilotcontamination_start;
X = [[xK;xM] [yK;yM]];
numcul = M/tau_cf;
%rng('default') % For reproducibility
opts = statset('Display','final');
[id,C] = kmeans(X,numcul,'Distance','cityblock',...
    'Replicates',numcul,'Options',opts);
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
       if(a<=tau_cf)
          matrix =  randperm(tau_cf,a);
       else
          m =  randperm(tau_cf,tau_cf);
          matrix = [m randperm(tau_cf,a-tau_cf)];
       end
       for k=1:a
           idx1(l(k)) = matrix(k);
       end
    end
end
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

%random pilot assignment

% rateavr_random = mean(Rate_start); pc_random = Pilotcontamination_start;
% for i=1:iteration
% randompilot = randi([1 tau_cf],1,K); %random greedy
% PilotSet_random = [];
% for j=1:K
%        %PilotSet caculates rate
%        PilotSet_random = [PilotSet_random pilot(:,randompilot(j))];
% end
% [C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet_random, nbrOfRealizations);
% Rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet_random, nbrOfRealizations);
% rateavr_random = [rateavr_random mean(Rate)];
% pc_random = [pc_random pilot_contamination(PilotSet_random,Beta,K,M,nbrOfRealizations)];
% end
% 
% %vet can minpilotcontamination
% matrix = [];
% for i=1:K
%    repeat = 4^(i-1);
%    repeat_1 = repmat(1,repeat,1);
%    repeat_2 = repmat(2,repeat,1);
%    repeat_3 = repmat(3,repeat,1);
%    repeat_4 = repmat(4,repeat,1);
%    repeat_col = [repeat_1;repeat_2;repeat_3;repeat_4];
%    repeat_col = repmat(repeat_col,(4^K / length(repeat_col)),1);
%    matrix = [repeat_col,matrix];
% end
% fitness = [];
% for i=1:tau_cf^K
%   fit = fitness_contamination(matrix(i,:),pilot,Beta,K,M,nbrOfRealizations);
%   fitness = [fitness fit];
% end
% index_random = find(fitness == min(fitness));
% PilotSet_vetcan_pc = matrix(index_random(1),:);
% % vet can rate
% fitness = [];
% for i=1:tau_cf^K
%   fit =  fitness_rate(matrix(i,:),pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
%   fitness = [fitness fit];
% end
% index_random = find(fitness == max(fitness));
% PilotSet_vetcan_mr = matrix(index_random(1),:);


figure()
x = 0:(iteration);
%plot(x,Rate_genetic_pc)
hold on
plot(x,Rate_genetic_rate)
hold on
plot(x,Rate_Greedy_Avrg)
hold on
plot(x,Rate_Greedy_Avrg1)
%plot(x,rateavr_random)
hold on
%plot(x,max(fitness)*ones(1,iteration+1))
hold on
title(strcat(int2str(K),' Users',int2str(M),' APs Average Downlink Rate'))
xlabel('Number of iterations')
ylabel('Average Downlink Rate (Mbits/s)')
legend('genetic (max average rate)','Rate greedy','Rate Kmean-Cluster')
figure()
plot(x,rate_pilotcontam_genetic)
hold on
plot(x,pilotcontamination_greedy)
hold on
plot(x,pilotcontamination_greedy1)
%plot(x,pc_random)
hold on
title('Pilot Contamination ')
xlabel('Number of iterations')
ylabel('Pilot Contamination ')
legend('genetic (max average rate)','Greedy','Kmean-Cluster')
%fitness_contamination(test,pilot,Beta,K,M,nbrOfRealizations)


