function [SC_GPASet PilotSet] = functionSC_GPA(M,K,RateEq42_dK,Krandomorder,mK_AP,Beta,PilotSet,nbrOfRealizations)
% Prepare for store Small Cell Greedy Pilot Assignment Set
SC_GPASet = zeros(M,K,nbrOfRealizations);
% Find the user with the lowest rate
[val kmin] = min(RateEq42_dK);
% Calculate the variance of the pilot contamination effect Algorithm 1
for nConta = 1:nbrOfRealizations
    subKConta = zeros(K,nbrOfRealizations);
    PilotContamination = zeros(K,nbrOfRealizations);
    for kchoosen = 1:K
            for kcomma = 1:K
                if (kcomma == kchoosen)
                    subKConta(kcomma,nConta) = 0;
                else
                    mkcomma = mK_AP(kcomma); % extract mkcomma AP serves kcomma-th user from mK_AP vector
                    subKConta(kcomma,nConta) = Beta(mkcomma,kchoosen,nConta) * abs(PilotSet(:,kchoosen,nConta)'...
                        * PilotSet(:,kcomma,nConta)).^2;
                end
            end
        PilotContamination(kchoosen) = sum(subKConta(:,nConta));
    end
    [val2 Pilotnew] = min(PilotContamination);
    PilotSet(:,kmin,nConta) = PilotSet(:,Pilotnew,nConta);
    SC_GPASet = PilotSet;
end
end