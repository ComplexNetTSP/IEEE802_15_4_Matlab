function NoiseFloor = CalNoiseFloor(B, NoiseFigure)
%% CalNoiseFloor 
% B : bandwidth in hz
% NoiseFigure : Noise Figure in dB
% NoiseFloor in dBm

% Boltzman Constant
BoltzmanCts = 1.3806504*10^(-23);
% Temperature in degree K (for 20 Degree Celcius)
T = 300;
% Noise Floor in dBm
NoiseFloor = 10*log10(BoltzmanCts * T * B * 10^3) + NoiseFigure + 1; 