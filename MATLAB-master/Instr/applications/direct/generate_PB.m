function generate_PB(PB, sampleRateDAC, inst)
%need to have them divisible by 32
% keys(PB): indicates channel
% values(PB): 2D array where the first column indicate 0s or 1s and 2nd
% column indicates the duration of those pulses.
    function setTask_PB(ch, num_ms_l, on_or_off, is_mult_ms, inst)
        inst.SendScpi(sprintf(':INST:CHAN %d',ch));
        inst.SendScpi('TASK:ZERO:ALL');
        %number of tasks: 2 from initial and final holding segment
        num_tasks = length(find(num_ms_l))+length(num_ms_l) + 2;
        %task index: t_idx
        t_idx = 1;
        % inital holding segment
        inst.SendScpi(sprintf(':TASK:COMP:LENG %d', num_tasks));
        inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
        t_idx = t_idx + 1;
        inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
        inst.SendScpi(':TASK:COMP:ENAB CPU');
        inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',3));
        inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
        inst.SendScpi(':TASK:COMP:TYPE SING');
        seg_idx = 1;
        for idx = 1:length(num_ms_l)
            if num_ms_l(idx) > 0
                inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
                t_idx = t_idx + 1;
                if on_or_off(idx) == 0
                    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', 1));
                    fprintf("assign task with repetitive positive 1ms \n");
                else
                    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', 2));
                    fprintf("assign task with repetitive negative 1ms \n");
                end
                inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', num_ms_l(idx)));
                inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
                inst.SendScpi(':TASK:COMP:TYPE SING');
            end
            %whenever there is a need for an additional segment
            if is_mult_ms(idx) == 0
                inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
                t_idx = t_idx + 1;
                inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', seg_idx+3));
                inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', 1));
                inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
                inst.SendScpi(':TASK:COMP:TYPE SING');
                seg_idx = seg_idx + 1; 
            end
        end
        inst.SendScpi(sprintf(':TASK:COMP:SEL %d', t_idx));
        inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
        inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',3));
        inst.SendScpi(':TASK:COMP:TYPE SING');
        inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',1));
        inst.SendScpi('TASK:COMP:WRITE');
        resp = inst.SendScpi('SOUR:FUNC:MODE TASK');
    end
k = keys(PB);
assert(length(PB) == 1, "Tabor PB only supports one channel");
for i = 1: length(PB)
    %get channel and marker number
    ch = k{i};
    PSeq = PB(ch);
    % num_ms_l indicates how much 1ms each low or high signal contains
    num_ms_l = [];
    % on_or_off indicates whether the signal from the marker is low or high
    on_or_off = [];
    % is_mult_ms indicates whether the low or high signal is multiple of
    % 1ms, if it isn't then download one more segment.
    is_mult_ms = [];
    %first store 1ms segment ON and OFF semgnet for the 1st marker
    segMem = 1;
    [I, Q] = makeDC(1e-3*sampleRateDAC);
    seg_1ms_off = uint8(zeros(1, length(I)));
    seg_1ms_on = uint8(ones(1, length(I)));
    %create segment with 1ms off
    downLoadIQ(ch, segMem, I, Q, inst);
    downLoad_mrkr(ch, segMem, seg_1ms_off, seg_1ms_off, inst);
    segMem = segMem + 1;
    %create segment with 1ms on
    downLoadIQ(ch, segMem, I, Q, inst);
    downLoad_mrkr(ch, segMem, seg_1ms_on, seg_1ms_off, inst);
    segMem = segMem + 1;
    %create markhold for initial and final pulses
    DClen = 64;
    [holdI, holdQ] = makeDC(DClen);
    markHold = uint8(zeros(1, DClen));
    downLoadIQ(ch, segMem, holdI, holdQ, inst);
    downLoad_mrkr(ch, segMem, markHold, markHold, inst);
    segMem = segMem + 1;
    for j = 1 : length(PSeq)
        on_or_off(end+1) = PSeq(j,1);
        mrker_time = PSeq(j,2);
        if mrker_time >= 1e-3
            num_ms = floor(mrker_time/1e-3);
            num_ms_l(end+1) = num_ms;
            mrker_time = mrker_time - num_ms*1e-3;
        else
            num_ms_l(end+1) = 0;
        end
        if mrker_time == 0
            is_mult_ms(end+1) = 1;
        else
            is_mult_ms(end+1) = 0;
            len_seg = 64*round(sampleRateDAC*mrker_time/64);
            new_seg = uint8(PSeq(j,1)) * uint8(ones(1, len_seg));
            mrkr2_seg = uint8(zeros(1, len_seg));
            [I, Q] = makeDC(len_seg);
            downLoadIQ(ch, segMem, I, Q, inst);
            downLoad_mrkr(ch, segMem, new_seg, mrkr2_seg, inst);
            segMem = segMem + 1;
        end
    end
    setTask_PB(ch, num_ms_l, on_or_off, is_mult_ms, inst);
end

end