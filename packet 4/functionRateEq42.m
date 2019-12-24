function RateEq42_dK = functionRateEq42(M, K, Gchannel, Krandomorder, mK_AP, Beta, rhod_sc, rhodp_sc, taud_sc, PilotSet, nbrOfRealizations)
% function [Alpha_dK Nuy_MK MeanNuy_MK RateEq42_dK] = functionRateEq42(M, K, Gchannel, Krandomorder, mK_AP, Beta, rhod_sc, rhodp_sc, taud_sc, PilotSet, nbrOfRealizations)


% INPUT
% We assume that each user is served by only one AP. Therefore, K-APs
% available denotes by mK_AP matrix will serve K-users in Krandomorder
% matrix corresponding
% M                 = number of APs
% K                 = number of users
% Gchannel          = M x K x nbrOfRealizations
% Krandomorder      = 1 x K x nbrOfRealizations
% mK_AP             = 1 x K x nbrOfRealizations
% Beta              = M x K x nbrOfRealizations Large Scale Fading matrix
% rhod_sc           = Normalized downlink transmit power of Small Cell
% rhodp_sc          = Normalizaed downlink pilot transmit power Small Cell
% taud_sc           = length of downlink pilot sequence in Small Cell
% PilotSet          = tau_cf x K x nbrOfRealizations Pilot matrix (tau_cf =
% taud_sc = tauu_sc)
% nbrOfRealizations 

% OUTPUT
% Alpha_dK          = 1 x K vector power control coefficient
% Nuy_MK            = K x K x nbrOfRealizations
% MeanNuy_MK        = K x K x nbrOfRealizations
% RateEq42          = K x nbrOfRealizations
Alpha_dK = ones(1,K);
% Calculate the covariance matrix of the MMSE estimate of g_mk channel
% Nuy_MK 
Nuy_MK = zeros(K,K,nbrOfRealizations);
Numerator = zeros(K,K,nbrOfRealizations);
Denominator = zeros(K,K,nbrOfRealizations);

MeanNuy_MK = zeros(K,K,nbrOfRealizations);
Mean_Numerator = zeros(K,K,nbrOfRealizations);
Mean_Denominator_2nd = zeros(K,K,nbrOfRealizations);
Mean_Denominator = zeros(K,K,nbrOfRealizations);

Numerator42 = zeros(K,K,nbrOfRealizations);
RateEquation42 = zeros(K,K,nbrOfRealizations);
for nSCRate = 1:nbrOfRealizations
    for kSCRate = 1:length(Krandomorder)
        kuser = Krandomorder(kSCRate);
        for mSCRate = 1:length(mK_AP)
            subDenominator = zeros(1,K); % Calculate Eq(40)
            subMean_Denominator = zeros(1,K); % Calculate Eq (44)
            
            % Calculate Downlink Rate for user k-th, which is served by
            % mkk-th AP. Note that: matrix Krandomorder and mK_AP store
            % a number of user/AP but no order.
            if (mSCRate == kSCRate)
                mkAP = mK_AP(mSCRate);
                Numerator(mSCRate,kSCRate,nSCRate) = taud_sc * rhodp_sc * Beta(mkAP,kuser,nSCRate).^2;
                for kcomma = 1:K
                    kindex = find(Krandomorder == kcomma);
                    m_kcomma = mK_AP(kindex);
                    subDenominator(1,kcomma) = Beta(m_kcomma,kuser,nSCRate) * ...
                        (abs(PilotSet(:,kuser,nSCRate)' * PilotSet(:,kcomma,nSCRate))).^2;
                end
                Denominator(mSCRate,kSCRate,nSCRate) = taud_sc * rhodp_sc * sum(subDenominator,2) + 1;
                Nuy_MK(mSCRate,kSCRate,nSCRate) = Numerator(mSCRate,kSCRate,nSCRate)./Denominator(mSCRate,kSCRate,nSCRate); %Eq (40)
                Mean_Numerator(mSCRate,kSCRate,nSCRate) = rhod_sc * Alpha_dK(kuser) * Nuy_MK(mSCRate,kSCRate,nSCRate);
                for kkcomma = 1:K
                    kkindex = find(Krandomorder == kkcomma);
                    m_kkcomma = mK_AP(kkindex);
                    if (kkcomma ~= kuser)
                        subMean_Denominator(1,kkcomma) = Alpha_dK(kkcomma) * Beta(m_kkcomma,kuser,nSCRate);
                    else 
                        subMean_Denominator(1,kkcomma) = 0;
                    end
                end
                Mean_Denominator_2nd(mSCRate,kSCRate,nSCRate) = rhod_sc * sum(subMean_Denominator,2);
                Mean_Denominator(mSCRate,kSCRate,nSCRate) = 1 + rhod_sc * Alpha_dK(kuser) *...
                    (Beta(mkAP,kuser,nSCRate) - Nuy_MK(mSCRate,kSCRate,nSCRate))...
                    + Mean_Denominator_2nd(mSCRate,kSCRate,nSCRate);
                Numerator42(mSCRate,kSCRate,nSCRate) =  rhod_sc * Alpha_dK(kuser) * abs(Gchannel(mkAP,kuser,nSCRate)).^2;
%                 MeanNuy_MK(mSCRate,kSCRate,nSCRate) = Mean_Numerator(mSCRate,kSCRate,nSCRate)./Mean_Denominator(mSCRate,kSCRate,nSCRate); % Eq (44)
                % Calculate Downlink Rate by Eq (42). We use perfectly
                % effective channel gain (perfect channel state information
                % was estimated by uplink transmission
                RateEquation42(mSCRate,kSCRate,nSCRate) = log2(1 + Numerator42(mSCRate,kSCRate,nSCRate)./...
                    Mean_Denominator(mSCRate,kSCRate,nSCRate));
            end
        end
    end
    RateEq42_dK(:,nSCRate) = diag(RateEquation42);
end
end