function EbNo = RSSIToEbNo(RSSI,NoiseFloor,DataRate,Bw,codeRate)
%% RSSIToEbNo
% SNR : in dB
% RSSI : in dBm
% Noise floor : in dBm
% DataRate : in bits
% nBits : Number of bits per symbol
% Bw : Bandwidth in hz
% codeRate : WIFI (1/11)

%EbNo = ( 10.^((RSSI - NoiseFloor)./10) ) ./ ((DataRate*codeRate)/Bw);
Gain = 10*log10((DataRate*codeRate)/Bw);
EbNo_dBm = RSSI - NoiseFloor - Gain;
%fprintf('ebnodbm = %f\n',EbNo_dBm);
EbNo = dBm2mW(EbNo_dBm);
%fprintf('ebno_temp = %f\n',EbNo);