function [rate_avrg] = fitness_rate(a,pilot,M, K, PowerControl,Beta,tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations)
PilotSet = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet = [PilotSet pilot(:,a(j))];
end
% code rate
[PilotContamination_avgr, subMConta1] = fitness_contamination(a,pilot,Beta,K,M,nbrOfRealizations);
[C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
RateEq24_dK = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
rate_avrg = mean(RateEq24_dK);










