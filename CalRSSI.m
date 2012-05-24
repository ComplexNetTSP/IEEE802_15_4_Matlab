function rssi = CalRSSI_RIMmodel(d, alpha, lambda, ShadowingStd, Prdbm, Noise, d0, distmin, distmax)
% alpha Attenuation Factor
%   alpha = 4 => Indoor
%   alpha = 2 => Free Space
%
% Prdbm : Transmission power in dBm 
%   puissance d'emmision en dBm 100mW => 20dBm
%   Puissance en dBm = 10*log10(Puissance en wath) + 10 * log10(10^3)
%   example :  
%   Prdbm = 10*log10(100*10^(-3)) + 10 * log10(10^3)
%
% Lambda pour le wifi 802.11b 2.4Ghz => 12.5 cm
%
% Fading en db
%%
% Pour le WiFi
PL_d0 = 20*log10((4*pi)/lambda) + Noise;
L = (PL_d0 + 10*alpha.*log10(d/d0))+ ShadowingStd;
%fprintf('>>>>>> L = %f\n',L);
rssi =  (Prdbm - L);