function [Hchannel Gchannel Wnoise] = functionGchannelGenerating(M, K, tau_cf, Beta, nbrOfRealizations)

% Store M x K x nbrOfRealizations channel matrix
Gchannel = zeros(M, K, nbrOfRealizations);
%Wnoise = zeros(tau_cf,M,nbrOfRealizations);

% Creat M x K x nbrOfRealizations small scale fading channel matrix
Hchannel = sqrt(0.5)*(randn(M, K, nbrOfRealizations) + 1i * randn(M, K, nbrOfRealizations));

% Creat tau_cf x M x nbrOfRealizations vector noise
Wnoise = sqrt(0.5)*(randn(tau_cf, M, nbrOfRealizations) + 1i * randn(tau_cf, M, nbrOfRealizations));
for nG = 1:nbrOfRealizations
    for mG = 1:M
        for kG = 1:K
            Gchannel(mG, kG, nG) = sqrt(Beta(mG, kG, nG)) * Hchannel(mG, kG, nG);
        end
    end
end
