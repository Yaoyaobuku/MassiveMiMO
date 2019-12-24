function [rate_avrg] = fitness_rate(a,pilot)
tau_cf= 4;
K=10; M=20; nbrOfRealizations = 1; D_sqr = 1000;
taud_sc = 20; tauu_sc = 20; BW = 20e6; NF_dB = 9;
AVErhod_cf = 200; AVErhou_cf = 100; AVErhop_cf = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DistanceControl = 'Uni'; % Control two solutions creat uniformly distributed
% 'Halton' use Halton sequenc, and 'Uni' use makdedist Uniformly-distribution
ShadowingControl = 'uncorrelated'; % Control two shadowing correlation model: 'uncorrelated' or 'correlated'
PowerControl = 'No'; % Two Power Control Mode: 'No' = without Power Control / 'Yes' = Max-Min Power Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d_MK xM yM xK yK] = functionDistance(M, K, D_sqr, DistanceControl, nbrOfRealizations);
[Beta PL z_MK] = functionLargeScaleFading(d_MK, M, K, ShadowingControl, nbrOfRealizations);
 % Beta = ones(M, K, nbrOfRealizations); % beta_mk = 1
[NoisePower rhod_cf rhou_cf rhop_cf rhod_sc rhou_sc rhoup_sc rhodp_sc] = functionNormalizedTransmitSNRs(M, K, BW, NF_dB, AVErhod_cf, AVErhou_cf, AVErhop_cf);
[Hchannel Gchannel Wnoise] = functionGchannelGenerating(M, K, tau_cf, Beta, nbrOfRealizations);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pilot = [[1 0 0 0];[0 1 0 0 ];[0 0 1 0];[0 0 0 0]] %generate pilot
PilotSet = [];
for j=1:10
       %PilotSet caculates rate
       PilotSet = [PilotSet pilot(:,a(j))];
end
   % code rate
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
RateEq24_dK = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
rate_avrg = min(RateEq24_dk);
[GreedyPilotSet PilotSet] = functionGreedyPilotAssignment(M,K,RateEq24_dK,Beta,PilotSet,nbrOfRealizations)
Rate_greedy = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, GreedyPilotSet, nbrOfRealizations);













