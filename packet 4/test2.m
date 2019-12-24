pilot_genetic_1 = [1 2 3 4 2 3 4 1 1 5 6 2];
PilotSet1 = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet1 = [PilotSet1 pilot(:,pilot_genetic_1(j))];
end
%avr_pilotcontam_genetic = [avr_pilotcontam_genetic pilot_contamination(PilotSet1,Beta,K,M,nbrOfRealizations)];
[C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet1, nbrOfRealizations);
Rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet1, nbrOfRealizations);
mean(Rate)