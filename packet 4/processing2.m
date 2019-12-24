M=20; K=8; nbrOfRealizations = 1; D_sqr = 1000;
tau_cf = 4; taud_sc = 4; tauu_sc = 4; BW = 20e6; NF_dB = 9;
AVErhod_cf = 200; AVErhou_cf = 100; AVErhop_cf = 100;tau_c =200;
% nbrOfRealizations is loop to calculate mean in Eq (26) and Eq (42) when
% the user knows only the channel statistics. Therefore nbrOfRealizations
% is very large (2e6 - 4e6) but this paper analyzed to approximate Eq (42)
% by Eq (43). In this version code, we consider that the user perfectly
% knows its effective channel gain and we set nbrOfRealizations = 1.


%% Greedy Pilot Assignment. Note that: Set Nstep == 1 to run Random Pilot Assignment Mode 
Nstep = 7; % Nstep = the number of iterations K/4

NumMonteCarlo = 100; %200
DistanceControl = 'Uni'; % Control twdso solutions creat uniformly distributed
% 'Halton' use Halton sequenc, and 'Uni' use makdedist Uniformly-distribution
ShadowingControl = 'uncorrelated'; % Control two shadowing correlation model: 'uncorrelated' or 'correlated'
PowerControl = 'No'; % Two Power Control Mode: 'No' = without Power Control / 'Yes' = Max-Min Power Control
Rate_CF = zeros(K,NumMonteCarlo);
Rate_SC = zeros(K,NumMonteCarlo);
for Nloop = 1:NumMonteCarlo
    Nloop
    ncount = 1; % counter of Cell Free Greedy Pilot Assignment
    sc_ncount = 1; % counter of Small Cell Greedy Pilot Assignment
    RateEq24 = zeros(K,Nstep);
    RateEq42 = zeros(K,Nstep);
    
    [d_MK xM yM xK yK] = functionDistance(M, K, D_sqr, DistanceControl, nbrOfRealizations);
    [Beta PL z_MK] = functionLargeScaleFading(d_MK, M, K, ShadowingControl, nbrOfRealizations);
    % Beta = ones(M, K, nbrOfRealizations); % beta_mk = 1
    [NoisePower rhod_cf rhou_cf rhop_cf rhod_sc rhou_sc rhoup_sc rhodp_sc] = functionNormalizedTransmitSNRs(M, K, BW, NF_dB, AVErhod_cf, AVErhou_cf, AVErhop_cf);
    [Hchannel Gchannel Wnoise] = functionGchannelGenerating(M, K, tau_cf, Beta, nbrOfRealizations);
    PilotSet = functionRandomPilotAssignment(tau_cf, K, nbrOfRealizations);
    %% SMALL CELL
    [Krandomorder mK_AP] = functionAPSelection(M,K,Beta,nbrOfRealizations);
    while (sc_ncount <= Nstep)
        RateEq42_dK = functionRateEq42(M, K, Gchannel, Krandomorder, mK_AP, Beta, rhod_sc, rhodp_sc, taud_sc, PilotSet, nbrOfRealizations);
        RateEq42(:,sc_ncount) = RateEq42_dK;
        [SC_GPASet PilotSet] = functionSC_GPA(M,K,RateEq42_dK,Krandomorder,mK_AP,Beta,PilotSet,nbrOfRealizations);
        sc_ncount = sc_ncount + 1;
    end
    Rate_SC(:,Nloop) = RateEq42_dK;
    %% CELL-FREE
    while (ncount <= Nstep)
        % [C, Cw, Gest, Gamma, GammaMean, Eta, EtaMean] = functionPerfectMMSE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Beta, Wnoise, PilotSet, nbrOfRealizations);
        % % [C, Gest, Gamma, GammaMean Eta EtaMean] = functionMMSECE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, Beta, PilotSet, nbrOfRealizations);
        % % [C, Gest, Gamma, Eta] = functionMMSE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Beta, Wnoise, PilotSet, nbrOfRealizations);
        [C, Gest, Gamma, Eta] = functionCE(M, K, PowerControl, tau_cf, rhop_cf, Gchannel, Wnoise, PilotSet, nbrOfRealizations);
        RateEq24_dK = functionCalculateRateEq24(M, K, rhod_cf, Eta, Gamma, Beta, PilotSet, nbrOfRealizations);
        RateEq24(:,ncount) = RateEq24_dK;
        [GreedyPilotSet PilotSet] = functionGreedyPilotAssignment(M,K,RateEq24_dK,Beta,PilotSet,nbrOfRealizations);
        ncount = ncount + 1;
    end
    Rate_CF(:,Nloop) = RateEq24_dK;
    clear d_MK xM yM xK yK Beta PL z_MK Hchannel Gchannel Wnoise PiloSet C Gest Gamma Eta RateEq24_dK RateEq24 RateEq42 Krandomorder mK_AP
end

