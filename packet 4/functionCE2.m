function [C, Gest, Gamma, Eta] = functionCE2(M, K,Beta,subMConta1, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations)
% INPUT
% M             = Number of APs
% K             = Number of users
% tau_cf        = Length of pilot sequence
% rhop_cf       = Normalized transmit power of pilot symbols
% Gchannel      = M x K x nbrOfRealizations channel matrix
% Wnoise        = tau_cf x M x nbrOfRealizations vector noise
% PilotSet      = tau_cf x K x nbrOfRealizations random pilot matrix
% Iter          = the loop to calculate Mean in (4) Equation
% nbrOfRealizations = Number of channel realizations
% OUTPUT
% C             = M x K value in equation (4)
% Gest          = M x K x nbrOfRealizations estimated channel matrix
% Gamma; Eta    = M x K matrix in equation (8)

% Prepare to store OUTPUT
C = zeros(M,K,nbrOfRealizations);
Gest = zeros(M,K,nbrOfRealizations);
Gamma = zeros(M,K,nbrOfRealizations);

for nCE = 1:nbrOfRealizations
    
    % Prepare to store tau_cf x 1 received pilot vector at M APs. Therefore
    % y_pm = (tau_cf x M) matrix
    y_pm = zeros(tau_cf,M);
    for mCE = 1:M
        
        % Prepare to store y_pm per K pilot sequence
        y_pmsub = zeros(tau_cf,K);

        for kCE = 1:K
            y_pmsub(:,kCE) = Gchannel(mCE,kCE,nCE) * PilotSet(:,kCE,nCE);
            if (kCE == K) 
                y_pm(:,mCE) = sqrt(tau_cf*rhop_cf)*sum(y_pmsub,2) + Wnoise(:,mCE,nCE); %Eq (2)
            end
        end
    end
    Yest_pmk = zeros(M, K, nbrOfRealizations);
    G_numeratormk = zeros(M, K, nbrOfRealizations);
    Yest_abs = zeros(M, K, nbrOfRealizations);
    for mm = 1:M
        for kk = 1:K
            Yest_pmk(mm,kk,nCE) = PilotSet(:,kk)' * y_pm(:,mm); %Eq (3)
            C(mm,kk,nCE) = (sqrt(tau_cf*rhop_cf)*Beta(mm,kk))/(tau_cf*rhop_cf*(subMConta1(mm,kk)+Beta(mm,kk)) +1);
        end
    end
end
%Calculate C: M x K and Gest: M x K x nbrOfRealizations
for nC = 1:nbrOfRealizations
    for mC = 1:M
        for kC = 1:K
            Gest(mC,kC,nC) = C(mC,kC)*Yest_pmk(mC,kC,nC); %Eq (4)
        end
    end
end

%Calculate Gamma: M x K and Eta: M x K
for mGamma = 1: M
    for kGamma = 1:K
        Gamma(mGamma,kGamma) =(1/nbrOfRealizations)*...
            norm(reshape(Gest(mGamma,kGamma,:),[nbrOfRealizations,1])).^2; % Equation (8)
    end
end
if (strcmp(PowerControl, 'No') == 1)
    for mEta = 1:M
        Eta(mEta,:) = (1/sum(Gamma(mEta,:))) * ones(1, K);
    end
end