function rate = fitnessRate_GA(a, pilot)
PilotSet = [];
for j=1:K
     %PilotSet caculates rate
     PilotSet = [PilotSet pilot(:,a(j))];
end
   % code rate
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);