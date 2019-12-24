function [PilotContamination_avgr] = pilot_contamination(pilot,Beta,K,M,nbrOfRealizations)
PilotSet = pilot;
for nConta = 1:nbrOfRealizations
subMConta = zeros(M,nbrOfRealizations);
subKConta = zeros(M,K,nbrOfRealizations);
PilotContamination = zeros(K,nbrOfRealizations);
for kchoosen = 1:K
    for mConta = 1:M
        for kConta = 1:K
            if (kConta == kchoosen)
                    subKConta(mConta,kConta,nConta) = 0;
            else
                    subKConta(mConta,kConta,nConta) = Beta(mConta,kConta,nConta) * abs(PilotSet(:,kchoosen,nConta)'...
                        * PilotSet(:,kConta,nConta)).^2;
            end
        end
            subMConta(mConta,nConta) = sum(subKConta(mConta,:,nConta),2);
    end
PilotContamination(kchoosen) = sum(subMConta(:,nConta));
end
end
PilotContamination_avgr  = mean(PilotContamination);