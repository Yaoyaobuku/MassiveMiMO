function [Beta PL z_MK] = functionLargeScaleFading(d_MK, M, K, ShadowingControl, nbrOfRealizations)
%%Large Scale Fading Model includes Path loss Model and Shadowing
%%Correlation Model. This code version employs the Hata-COST231 propagation
%%model. Therefore, d1 and d0 are the reference distance

% INPUT (TABLE I)
% d_MK: the (M x K x nbrOfRealizations) matrix distance between m-th AP and k-th user (m)
% FreqCarrier (MHz): the carrier frequency (MHz) = 1900 GHz
% h_AP (m): the AP antenna height = 15 m 
% h_u (m): the user antenna height = 1.65 m
% d0 < d1 (m): the reference distances (10m & 50m)
% rho_shadow (dB): the standard deviation of the shadowing = 8 dB
% L (dB): in the (53) equation, L is a constant piece of path loss model, which
% depends on the carrier frequency, the AP antenna height, the user antenna
% height.

% OUTPUT
% Beta (absolute): M x K x nbrOfRealizations Large Scale Fading matrix
% PL (dB): M x K x nbrOfRealizations path loss model matrix
% z_MK (absolute): M x K x nbrOfRealizations Shadowing Correlation
% Coefficients matrix

FreqCarrier = 1900;
h_AP = 15; 
h_u = 1.65;
d0 = 10; 
d1 = 50;
L = 46.3 + 33.9*log10(FreqCarrier)- 13.82*log10(h_AP)-(1.1*log10(FreqCarrier)-0.7)*h_u+(1.56*log10(FreqCarrier)-0.8);
rho_shadow = 8;

% Creat the (M x K x nbrOfRealizations) matrix for storing the large scale
% fading
Beta = zeros(M, K, nbrOfRealizations);

% Creat the (M x K x nbrOfRealizations) matrix for storing the Path loss
% model
PL = zeros(M, K, nbrOfRealizations);

% Creat the (M x K x nbrOfRealizations) matrix for storing the Shadowing
% Coefficients
z_MK = zeros(M, K, nbrOfRealizations);

for nLSF = 1:nbrOfRealizations
    
    % Calculate the Shadowing Correlation Model z_MK
    % Note that: There is no shadowing when d_MK(m, k , nbrOfRealizations) < d1
    if (strcmp(ShadowingControl, 'uncorrelated') == 1)
        for mShadow = 1:M
            for kShadow = 1:K
                if d_MK(mShadow, kShadow, nLSF) > d1
                    z_MK(mShadow, kShadow, nLSF) = randn;
                else d_MK(mShadow, kShadow, nLSF) <= d1
                    z_MK(mShadow, kShadow, nLSF) = 0;
                end
            end
        end
    end
    % Calculate the path loss model PL_mk (dB) by equation (52) and Large
    % Scale Fading (absolute value) by equation (51)
    % Note that: In equation (52), we have change d_mk unit to (km) 
    for mLSF = 1:M
        for kLSF = 1:K
            if d_MK(mLSF, kLSF, nLSF) > d1
                PL(mLSF, kLSF, nLSF) = -L - 35*log10(0.001 * d_MK(mLSF, kLSF, nLSF));
            elseif d0 < d_MK(mLSF, kLSF, nLSF) <= d1
                PL(mLSF, kLSF, nLSF) = -L - 15*log10(0.001 * d1) - ...
                    20*log10(0.001 * d_MK(mLSF, kLSF, nLSF));
            else d_MK(mLSF, kLSF, nLSF) <= d0
                PL(mLSF, kLSF, nLSF) = -L - 15*log10(0.001 * d1) - ...
                    20*log10(0.001 * d0);
            end
            Beta(mLSF, kLSF, nLSF) = 10.^((PL(mLSF, kLSF, nLSF) + ...
                rho_shadow * z_MK(mLSF, kLSF, nLSF))./10);
        end
    end
end