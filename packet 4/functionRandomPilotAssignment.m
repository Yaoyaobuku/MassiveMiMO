function PilotSet = functionRandomPilotAssignment(tau_cf, K, nbrOfRealizations)

%Prepare to store K Pilot sequences set
PilotSet = zeros(tau_cf, K, nbrOfRealizations);

%Prepare to generate Pilot sequences that has unit power vector
RandomSequence = sqrt(0.5) * (randn(tau_cf, K, nbrOfRealizations) + 1i*randn(tau_cf, K, nbrOfRealizations));
for nP = 1:nbrOfRealizations
    if (tau_cf >= K)
        [s v d] = svd(RandomSequence(:, :, nP));
        PilotSet(:,:,nP) = s(:,1:K);
    elseif (tau_cf < K)
        [S v d] = svd(RandomSequence(:, :, nP));
        PilotSet(:, 1:tau_cf, nP) = S;
        Xgain = K/tau_cf;
        count_tau_cf = tau_cf;
        while (Xgain > 2)
            PilotSet(:, (count_tau_cf + 1):(count_tau_cf + tau_cf), nP)...
                = S(:, randperm(tau_cf,tau_cf));
            Xgain = Xgain - 1;
            count_tau_cf = count_tau_cf + tau_cf;
        end
        if (Xgain <= 2)
            PilotSet(:, (count_tau_cf+1):K, nP) = S(:, randperm(tau_cf, K - count_tau_cf));
        end
    end
end
end