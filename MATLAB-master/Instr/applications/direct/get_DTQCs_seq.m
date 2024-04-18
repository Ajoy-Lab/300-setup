function [rem_l, seg_idx_l] = get_DTQCs_seq(sampleRateDAC, granularity, y_pulse_len, x_pulse_len, spacing, T1, T2, seq_time)
    % y-segment length
    y_seg_len = y_pulse_len + spacing;
    y_seg_len = round_to_DAC_freq(y_seg_len, sampleRateDAC, granularity);
    
    % x-segment length
    x_seg_len = x_pulse_len + spacing;
    x_seg_len = round_to_DAC_freq(x_seg_len, sampleRateDAC, granularity);
    
    time_arr = round_to_DAC_freq([T1, T2], sampleRateDAC, granularity);
    [T1, T2] = deal(time_arr(1), time_arr(2));
    [freq1, freq2] = deal(1/T1, 1/T2);
    Y1_l = [0, y_seg_len];
    Y2_l = [y_seg_len, 2*y_seg_len];
    
    % add all possible Y1 pulses
    while Y1_l(end, 2) + T1 < seq_time
        Y1_l = [Y1_l; Y1_l(end, 1)+T1, Y1_l(end, 2)+T1];
    end
    
    % add all possible Y2 pulses
    while Y2_l(end, 2) + T2 < seq_time
        Y2_l = [Y2_l; Y2_l(end, 1)+T2, Y2_l(end, 2)+T2];
    end
    
    % Check whether the two lists of y-pulses are incoherent
    Y_l = [Y1_l; Y2_l];
    Y_l = sortrows(Y_l, 1);
    overlap = false;
    for idx = (1: length(Y_l)-1)
        if Y_l(idx, 2) > Y_l(idx+1, 1)
            overlap = true;
            fprintf("%d, %d \n", Y_l(idx, 2), Y_l(idx+1, 1));
        end
    end
    assert(overlap == false, "Y-pulse lists need to be incoherent \n");
    
    % now get the number of x-pulse spacing and generate seg_idx_l to
    % generate the full pulse sequence
    num_x_pulse_l = [];
    % first segment is shortest DC segment
    [y_seg_idx, x_seg_idx] = deal(2, 3);
    % rem_seg_idx is segment index that corresponds to the last index
    rem_seg_idx = 4;
    
    seg_idx_l = [];
    rem_l = [];
    t_unit = granularity * vpa(1/sampleRateDAC);
    for idx = (1:length(Y_l) - 1)
        Y_spacing = Y_l(idx+1, 1) - Y_l(idx, 2);
        % add y-pulse
        seg_idx_l = [seg_idx_l; y_seg_idx, 1];
        % add x-pulse
        num = fix(Y_spacing/x_pulse_len);
        seg_idx_l = [seg_idx_l; x_seg_idx, num];
        
        % add remainder
        rem = mod(Y_spacing, x_seg_len);
        rem = round_to_DAC_freq(rem, sampleRateDAC, granularity);
        
        % if remainder is close to 0
        if rem >= t_unit
            seg_idx_l = [seg_idx_l; rem_seg_idx, 1];
            rem_l = [rem_l ; rem_seg_idx, rem];
            rem_seg_idx = rem_seg_idx + 1;
        end
    end
end