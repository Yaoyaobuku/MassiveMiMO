function RateEq24_dK = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations)
%PilotSet = pilot_genetic;
RateEq24_dK = zeros(K,nbrOfRealizations);
for nRate = 1:nbrOfRealizations
    Numerator = zeros(K, nbrOfRealizations);
    Denominator = zeros(K, nbrOfRealizations);
    for k0 = 1:K
        % Calculate Numerator of Eq (24)
        % Prepare to store subNum
        subNum = zeros(1,M);
        for mNum = 1:M
%             subNum(mNum) = sqrt(Eta(mNum,k0,nRate)) * Gamma(mNum,k0,nRate);
            subNum(mNum) = sqrt(Eta(mNum,k0)) * Gamma(mNum,k0);

        end
        Numerator(k0, nRate) = rhod_cf * (sum(subNum)).^2;
      
        subDeno_1st = zeros(1, K);
        subDeno_2nd = zeros(1, K);
        for kDeno = 1:K
         
            msubDeno_1st = zeros(1, M);
            msubDeno_2nd = zeros(1, M);
            if kDeno == k0
                subDeno_1st(kDeno) = 0; 
                subDeno_2nd(kDeno) = 0;
            else
                for mDeno = 1:M
%                     msubDeno_1st(mDeno) = sqrt(Eta(mDeno, kDeno, nRate)) * ...
%                         Gamma(mDeno, kDeno, nRate) * Beta(mDeno, k0, nRate) ./...
%                         Beta(mDeno, kDeno, nRate);
                    msubDeno_1st(mDeno) = sqrt(Eta(mDeno, kDeno)) * ...
                        Gamma(mDeno, kDeno) * Beta(mDeno, k0, nRate) ./...
                        Beta(mDeno, kDeno, nRate);
%                     msubDeno_2nd(mDeno) = Eta(mDeno, kDeno, nRate) * ...
%                         Gamma(mDeno, kDeno, nRate) * Beta(mDeno, k0, nRate);
                    msubDeno_2nd(mDeno) = Eta(mDeno, kDeno) * ...
                        Gamma(mDeno, kDeno) * Beta(mDeno, k0, nRate);
                end
                subDeno_1st(kDeno) = (sum(msubDeno_1st)).^2 * ((abs(PilotSet(:, kDeno, nRate)' *...
                    PilotSet(:, k0, nRate))).^2);
                subDeno_2nd(kDeno) = sum(msubDeno_2nd);
            end
        end
        Deno_1st = rhod_cf * sum(subDeno_1st);
        Deno_2nd = rhod_cf * sum(subDeno_2nd);
        Denominator(k0, nRate) = Deno_1st + Deno_2nd+1;
        RateEq24_dK(k0, nRate) = log2(1 + Numerator(k0, nRate)./Denominator(k0, nRate));
    end
end