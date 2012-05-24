%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2012 Telecom SudParis.
%
% Permission is hereby granted, free of charge, to any person obtaining
% a copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to
% permit persons to whom the Software is furnished to do so, subject to
% the following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB M-file CalRSSI_RIMmodel.m
% Authors : Mohamed-Haykel Zayani and Vincent Gauthier
% Emails : {mohamed-haykel.zayani, vincent.gauthier}@it-sudparis.eu
% Address : Laboratory CNRS S.A.M.O.V.A.R. - Dept RS2M
% Telecom Sud Paris * 9 rue C. Fourier * 91011 EVRY CEDEX * FRANCE
% Created : April 28th, 2010
% Updated : May 18th, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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