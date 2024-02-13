function generate_PB(PB, sampleRateDAC, inst)
%need to have them divisible by 32
% keys(PB): indicate channel and the number of a marker
% values(PB): 2D array where the first column indicate 0s or 1s and 2nd
% column indicates the duration of those pulses.
    function setTask_PB(ch, num_ms_l, on_or_off, inst)
        inst.SendScpi(sprintf(':INST:CHAN %d',ch));
        inst.SendScpi('TASK:ZERO:ALL');
        %number of tasks
        num_tasks = length(find(num_ms_l))+length(num_ms_l) + 2;
        %task index: t_idx
        t_idx = 1;
        inst.SendScpi(sprintf(':TASK:COMP:LENG %d', num_tasks)); % this should be more general?
        inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
        t_idx = t_idx + 1;
        inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
        inst.SendScpi(':TASK:COMP:ENAB CPU');
        inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',3));
        inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
        inst.SendScpi(':TASK:COMP:TYPE SING');
        for seg_idx = 1:length(num_ms_l)
            inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
            t_idx = t_idx + 1;
            if on_or_off(seg_idx) == 0
                inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', 1));
                fprintf("assign task with repetitive positive 1ms");
            else
                inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', 2));
                fprintf("assign task with repetitive negative 1ms");
            end
            inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', num_ms_l(seg_idx)));
            inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
            inst.SendScpi(':TASK:COMP:TYPE SING');
            inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
            t_idx = t_idx + 1;
            inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', seg_idx+3));
            inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', 1));
            inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
            inst.SendScpi(':TASK:COMP:TYPE SING');
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
    key = k{i};
    PSeq = PB(key);
    num_ms_l = [];
    on_or_off = [];
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
    for j = 1 : length(PB)
        mrker_time = PSeq(j,2);
        if mrker_time >= 1e-3
            num_ms = floor(mrker_time/1e-3);
            num_ms_l(end+1) = num_ms;
            mrker_time = mrker_time - num_ms*1e-3;
        else
            num_ms_l(end+1) = 0;
        end
        len_seg = 64*round(sampleRateDAC*mrker_time/64);
        new_seg = uint8(PSeq(j,1)) * uint8(ones(1, len_seg));
        on_or_off(end+1) = PSeq(j,1);
        mrkr2_seg = uint8(zeros(1, len_seg));
        [I, Q] = makeDC(len_seg);
        downLoadIQ(ch, segMem, I, Q, inst);
        downLoad_mrkr(ch, segMem, new_seg, mrkr2_seg, inst);
        segMem = segMem + 1;
    setTask_PB(ch, num_ms_l, on_or_off, inst);
    end
end

end