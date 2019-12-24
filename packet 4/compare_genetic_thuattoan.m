clear;
tau_cf= 10;
K=20; M=40; nbrOfRealizations = 1; D_sqr = 1000; population = 19;
taud_sc = 20; tauu_sc = 20; BW = 20e6; NF_dB = 9;
AVErhod_cf = 200; AVErhou_cf = 100; AVErhop_cf = 100; 
iteration = 30; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DistanceControl = 'Uni'; % Control two solutions creat uniformly distributed
% 'Halton' use Halton sequenc, and 'Uni' use makdedist Uniformly-distribution
ShadowingControl = 'uncorrelated'; % Control two shadowing correlation model: 'uncorrelated' or 'correlated'
PowerControl = 'No'; % Two Power Control Mode: 'No' = without Power Control / 'Yes' = Max-Min Power Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d_MK xM yM xK yK] = functionDistance(M, K, D_sqr, DistanceControl, nbrOfRealizations);
figure(1)
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
random = [1:tau_cf randi([1 tau_cf],1,K-tau_cf)]; %random greedy
PilotSet = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet = [PilotSet pilot(:,random(j))];
end
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
Rate_start = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
Pilotcontamination_start =  pilot_contamination(PilotSet,Beta,K,M,nbrOfRealizations);
%genentic pilot assigment with max average rate p1
pop = [];
PilotSet = [];
Rate_genetic_rate_p1 =mean(Rate_start);
rate_pilotcontam_genetic_p1 = Pilotcontamination_start;
pop = functioncreatpop(K,tau_cf,population,Beta);
pop = [random;pop];
for i=1:iteration
    parent = select_genetic(pop,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    child = cross_over_genetic(parent);
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
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet1, nbrOfRealizations);
Rate_genetic_rateavrg = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet1, nbrOfRealizations);
Rate_genetic_rate_p1 = [Rate_genetic_rate_p1 mean(Rate_genetic_rateavrg)];
rate_pilotcontam_genetic_p1 = [rate_pilotcontam_genetic_p1 pilot_contamination(PilotSet1,Beta,K,M,nbrOfRealizations)];
end

%genentic pilot assigment with max average rate p2
pop = [];
PilotSet = [];
Rate_genetic_rate_p2 =mean(Rate_start);
rate_pilotcontam_genetic_p2 = Pilotcontamination_start;
for i=1:(population)
    m = randi([1 tau_cf],1,K);
    pop = [pop;m];
end
pop = [random;pop];
for i=1:iteration
    parent = select_genetic(pop,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    child = cross_over_genetic(parent);
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
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet1, nbrOfRealizations);
Rate_genetic_rateavrg = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet1, nbrOfRealizations);
Rate_genetic_rate_p2 = [Rate_genetic_rate_p2 mean(Rate_genetic_rateavrg)];
rate_pilotcontam_genetic_p2 = [rate_pilotcontam_genetic_p2 pilot_contamination(PilotSet1,Beta,K,M,nbrOfRealizations)];
end

matrix = [];
for i=1:K
   % repeat = 4^(i-1);
   % repeat_1 = repmat(1,repeat,1);
 %   repeat_2 = repmat(2,repeat,1);
  %  repeat_3 = repmat(3,repeat,1);
  %  repeat_4 = repmat(4,repeat,1);
  %  repeat_col = [repeat_1;repeat_2;repeat_3;repeat_4];
  % repeat_col = repmat(repeat_col,(4^K / length(repeat_col)),1);
   % matrix = [repeat_col,matrix];
end
rate = [];
fitness = [];
%for i=1:1
 %  fit = fitness_rate(matrix(i,:),pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
  % fitness = [fitness fit];
%end
%max(fitness)

figure(2)
x = 0:(iteration);
plot(x,Rate_genetic_rate_p1)
hold on
plot(x,Rate_genetic_rate_p2)
%plot(x,max(fitness)*ones(1,iteration+1))
hold on
title('Average Downlink Rate')
xlabel('Number of iterations')
ylabel('Average Downlink Rate (Mbits/s)')
legend('genetic (max average rate) p1','genetic (max average rate) p2')
figure(3)
plot(x,rate_pilotcontam_genetic_p1)
hold on
plot(x,rate_pilotcontam_genetic_p2)
hold on
title('Pilot Contamination ')
xlabel('Number of iterations')
ylabel('Pilot Contamination ')
legend('genetic (max average rate) p1','genetic (max average rate) p2')

