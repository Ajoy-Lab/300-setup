function downLoadIQ(ch, segMem, dacWaveI, dacWaveQ, inst)
    fprintf(sprintf('Downloading waveform to channel %s, segment %s', num2str(ch), num2str(segMem)));

    dacWaveIQ = [dacWaveI; dacWaveQ];
    dacWaveIQ = dacWaveIQ(:)';
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi(':TRAC:FORM U16');
    inst.SendScpi(sprintf(':TRAC:DEF %d, %d',segMem, length(dacWaveIQ)));
    inst.SendScpi(sprintf(':TRAC:SEL %d',segMem));        
    % Download the binary data to segment
    prefix = ':TRAC:DATA 0,';
    %we must be using 16 bit system -- typecasting
    myWfm = uint16(dacWaveIQ);
    myWfm = typecast(myWfm, 'uint8');
    res = inst.WriteBinaryData(prefix, myWfm);
end