function downLoadIQ(ch, segMem, dacWaveI, dacWaveQ, inst)
    fprintf(sprintf('Downloading waveform to channel %s, segment %s \n', num2str(ch), num2str(segMem)));
%     disp('--- downLoadIQ called with arguments ---');
%     disp(['ch: ', num2str(ch)]);
%     disp(['segMem: ', num2str(segMem)]);
%     disp(['dacWaveI (length = ', num2str(length(dacWaveI)), '): ', mat2str(dacWaveI)]);
%     disp(['dacWaveQ (length = ', num2str(length(dacWaveQ)), '): ', mat2str(dacWaveQ)]);
%     disp(['inst: ', class(inst)]);

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