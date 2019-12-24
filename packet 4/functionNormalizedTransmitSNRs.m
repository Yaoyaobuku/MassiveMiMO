function [NoisePower rhod_cf rhou_cf rhop_cf rhod_sc rhou_sc rhoup_sc rhodp_sc] = functionNormalizedTransmitSNRs(M, K, BW, NF_dB, AVErhod_cf, AVErhou_cf, AVErhop_cf)

% INPUT
% M     : number of the APs
% K     : number of the users
% BW (Hz): 20 MHz
% NF_dB (dB): 9 dB Noise Figure
% AVErhod_cf (mW): 200 mW the transmit power of downlink data
% AVErhou_cf (mW): 100 mW the transmit power of uplink data
% AVErhop_cf (mW): 100 mW the transmit power of pilot symbols
% T0 = 290 (Kelvin): the noise temperature
% kB = 1.381*10^-23 (J/K): the Boltzmann constant

% OUTPUT
% NoisePower (W)
% rhod_cf (absolute): the normalized transmit SNR of downlink cell-free
% rhou_cf (absolute): the normalized transmit SNR of uplink cell-free
% rhop_cf (absolute): the normalized transmit SNR of pilot symbols
% rhod_sc (absolute): the normalizaed transmit SNR of downlink small-cell
% rhou_sc (absolute): the normalized transmit SNR of uplink small-cell

kB = 1.381*10^-23;
T0 = 290;
NoisePower = BW * kB * T0 * 10^(NF_dB/10);
rhod_cf = 0.001 * AVErhod_cf / NoisePower;
rhou_cf = 0.001 * AVErhou_cf / NoisePower;
rhop_cf = 0.001 * AVErhop_cf / NoisePower;
rhod_sc = M*rhod_cf/K;
rhou_sc = rhou_cf;
rhoup_sc = rhop_cf;
rhodp_sc = rhop_cf;
end