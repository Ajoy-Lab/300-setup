function downLoad_mrkr(ch, segMem, state1, state2, inst)
    fprintf('Downloading marker to channel %s, segment %s \n', num2str(ch), num2str(segMem))
%     disp(['ch: ', num2str(ch)]);
%     disp(['segMem: ', num2str(segMem)]);
%     disp(['state1 (length = ', num2str(length(state1)), '): ', mat2str(state1)]);
%     disp(['state2 (length = ', num2str(length(state2)), '): ', mat2str(state2)]);
    
    myMkr = uint8(state1 + 2*state2);
    
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi(':TRAC:FORM U16');
    inst.SendScpi(sprintf(':TRAC:SEL %d',segMem));
    
    myMkr = myMkr(1:2:length(myMkr)) + 16 * myMkr(2:2:length(myMkr));

    res = inst.WriteBinaryData(':MARK:DATA 0,', myMkr);
    assert(res.ErrCode == 0);
    
    inst.SendScpi(sprintf(':MARK:SEL %d',1));
    inst.SendScpi(':MARK:VOLT:PTOP 0.5');
    %inst.SendScpi(':MARK:VOLT:LEV 0.0')
    inst.SendScpi(':MARK:VOLT:OFFS 0.25');
    inst.SendScpi(':MARK:STAT ON');
    
    inst.SendScpi(sprintf(':MARK:SEL %d',2));
    inst.SendScpi(':MARK:VOLT:PTOP 1.0');
    %inst.SendScpi(':MARK:VOLT:LEV 0.0')
    inst.SendScpi(':MARK:VOLT:OFFS 0.0');
    inst.SendScpi(':MARK:STAT ON');
end