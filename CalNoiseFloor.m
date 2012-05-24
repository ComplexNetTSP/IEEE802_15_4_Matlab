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
% MATLAB M-file CalNoiseFloor.m
% Authors : Mohamed-Haykel Zayani and Vincent Gauthier
% Emails : {mohamed-haykel.zayani, vincent.gauthier}@it-sudparis.eu
% Address : Laboratory CNRS S.A.M.O.V.A.R. - Dept RS2M
% Telecom Sud Paris * 9 rue C. Fourier * 91011 EVRY CEDEX * FRANCE
% Created : April 28th, 2010
% Updated : May 18th, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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