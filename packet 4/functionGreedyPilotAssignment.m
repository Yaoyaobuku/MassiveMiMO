function [GreedyPilotSet PilotSet] = functionGreedyPilotAssignment(M,K,RateEq24_dK,Beta,pilot,PilotSet,nbrOfRealizations)
% Prepare for store Greedy Pilot Set
GreedyPilotSet = zeros(M,K,nbrOfRealizations);
l = size(pilot);
% Find the user with the lowest rate
[val kmin] = min(RateEq24_dK);
% Calculate the variance of the pilot contamination effect
for nConta = 1:nbrOfRealizations
    subMConta = zeros(M,nbrOfRealizations);
    subKConta = zeros(M,K,nbrOfRealizations);
    PilotContamination = zeros(l(1),nbrOfRealizations);
    for kchoosen = 1:l(1)
        for mConta = 1:M
            for kConta = 1:K
                if (kConta == kmin)
                    subKConta(mConta,kConta,nConta) = 0;
                else
                    subKConta(mConta,kConta,nConta) = Beta(mConta,kConta,nConta) * abs(pilot(:,kchoosen,nConta)'...
                        * PilotSet(:,kConta,nConta)).^2;
                end
            end
            subMConta(mConta,nConta) = sum(subKConta(mConta,:,nConta),2);
        end
        PilotContamination(kchoosen) = sum(subMConta(:,nConta));
    end
    [val2 Pilotnew] = min(PilotContamination);
    PilotSet(:,kmin,nConta) = pilot(:,Pilotnew,nConta);
    GreedyPilotSet = PilotSet;
end