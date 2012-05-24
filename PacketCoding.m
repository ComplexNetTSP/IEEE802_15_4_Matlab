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

function prrM = PacketCoding(pe,ENCODING,PREAMBLE_LENGTH,FRAME_LENGTH)
preseq = (1-pe).^(8*PREAMBLE_LENGTH);
if (ENCODING == 1)      % NRZ
    prrM = preseq.*((1-pe).^(8*(FRAME_LENGTH-PREAMBLE_LENGTH)));
elseif (ENCODING == 2)  % 4B5B
    prrM = preseq.*((1-pe).^(8*1.25*(FRAME_LENGTH-PREAMBLE_LENGTH)));
elseif (ENCODING == 3)  % MANCHESTER
    prrM = preseq.*((1-pe).^(8*2*(FRAME_LENGTH-PREAMBLE_LENGTH)));
elseif (ENCODING == 4)  % SECDED
    prrM = ((preseq.*((1-pe).^8)) + (8.*pe.*((1-pe).^7))).^((FRAME_LENGTH-PREAMBLE_LENGTH)*3);
else
    error('ENCODING is not correct');
end