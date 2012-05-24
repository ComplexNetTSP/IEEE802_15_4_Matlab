function pe = Perror(snr, MODULATION) 
if (MODULATION == 1)    % NCASK
    pe = 0.5*( exp(-0.5*snr) + Q( sqrt(snr) ) );
elseif(MODULATION == 2) % ASK
    pe = Q( sqrt(snr/2) );
elseif(MODULATION == 3) % NCFSK
    pe = 0.5*exp(-0.5*snr);
elseif(MODULATION == 4) % FSK
    pe = Q( sqrt(snr) );
elseif(MODULATION == 5) % BPSK
    pe = Q( sqrt(2*snr) );
elseif(MODULATION == 6) % DPSK
    pe = 0.5*exp(-snr);
else
    error('MODULATION is not correct');
end