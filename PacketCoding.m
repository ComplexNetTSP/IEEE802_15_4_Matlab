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