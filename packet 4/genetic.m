tau_cf= 6;
K=8; M=10; nbrOfRealizations = 1; D_sqr = 1000;
taud_sc = 20; tauu_sc = 20; BW = 20e6; NF_dB = 9;
AVErhod_cf = 200; AVErhou_cf = 100; AVErhop_cf = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DistanceControl = 'Uni'; % Control two solutions creat uniformly distributed
% 'Halton' use Halton sequenc, and 'Uni' use makdedist Uniformly-distribution
ShadowingControl = 'uncorrelated'; % Control two shadowing correlation model: 'uncorrelated' or 'correlated'
PowerControl = 'No'; % Two Power Control Mode: 'No' = without Power Control / 'Yes' = Max-Min Power Control
NumMonteCarlo = 1; %200
Rate_CF_genetic = zeros(K,NumMonteCarlo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Nloop = 1:NumMonteCarlo

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
PilotSet = functionRandomPilotAssignment(tau_cf, K, nbrOfRealizations);
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
Rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
pop = [];
PilotSet = [];
for i=1:10
    m = randi([1 tau_cf],1,K);
    pop = [pop;m];
end

for i=1:20
    parent = select_genetic(pop,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    child = cross_over_genetic(parent);
    pop = [parent;child];
end
fitness = [];
l = size(pop);
for i=1:l(1)
    fit = fitness_rate(pop(i,:),pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations);
    fitness = [fitness fit];
end
max(fitness)
index_genetic = find(fitness == max(fitness));
pilot_genetic = pop(index_genetic(1),:)
PilotSet = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet = [PilotSet pilot(:,pilot_genetic(j))];
end
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
Rate_genetic = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
Rate_CF_genetic(:,Nloop) = Rate_genetic;
mean(Rate_genetic);
clear d_MK xM yM xK yK Beta PL z_MK Hchannel Gchannel Wnoise PiloSet C Gest Gamma Eta RateEq24_dK RateEq24 RateEq42 Krandomorder mK_AP
end


