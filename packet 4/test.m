figure()
scatter(xM,yM,'^')
hold on
for i=1:length(xM) 
    text(xM(i)+10,yM(i)+10,int2str(i));
end
hold on
for i=1:length(xK) 
    scatter(xK(i),yK(i),'r');
    text(xK(i)+10,yK(i)+10,int2str(i));
end
PilotSet_vetcan_pc = pilot_genetic;

%fitness_rate(PilotSet_vetcan_pc,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations)
%fitness_rate(PilotSet_ventcan_mr,pilot,M, K, PowerControl,Beta, tau_cf,rhod_cf, rhop_cf, Gchannel, Wnoise, nbrOfRealizations)
PilotSet = [];
for j=1:K
       %PilotSet caculates rate
       PilotSet = [PilotSet pilot(:,PilotSet_vetcan_pc(j))];
end
[PilotContamination_avgr, subMConta1] = fitness_contamination_forgreedy(PilotSet,Beta,K,M,nbrOfRealizations);
[C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
rate = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
mean(rate)
figure()
scatter(xM,yM,'^');
hold on 
scatter(xK,yK,'r');
for i=1:length(xK) 
    text(xK(i)+10,yK(i)+10,strvcat(int2str(PilotSet_vetcan_pc(i)),strcat(num2str(rate(i)),'Mbits/s')));
end
title(strcat(int2str(K),' Users',int2str(M),' APs'));
xlabel('Distance 1 km2');
ylabel('Distance 1 km2');
legend('Access Point','User')













