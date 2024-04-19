function seg_idx_l = get_DTQCs_seq(sampleRateDAC, granularity, lengths, spacings, T1, T2, seq_time)
    % y-segment length
    fprintf("generating segment list for DTQC with period T1: %.6f, T2: %.6f ...\n",T1, T2);
    assert (lengths(2) == 2*lengths(1), "lengths(2) must be y-pulse");
    y_seg_len = lengths(2) + spacings(2);
    y_seg_len = round_to_DAC_freq(y_seg_len, sampleRateDAC, granularity);
    
    % x-segment length
    x_seg_len = lengths(3) + spacings(3);
    x_seg_len = round_to_DAC_freq(x_seg_len, sampleRateDAC, granularity);
    
    % x_seg_n_acq length
    x_seg_n_acq_len = lengths(4) + spacings(4);
    x_seg_n_acq_len = round_to_DAC_freq(x_seg_n_acq_len, sampleRateDAC, granularity);
    
    time_arr = round_to_DAC_freq([T1, T2], sampleRateDAC, granularity);
    [T1, T2] = deal(time_arr(1), time_arr(2));
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
    
    % first segment is initial (pi/2)_y pulse
    [y_seg_idx, x_seg_idx, x_seg_n_acq_idx, DC_unit_idx] = deal(2, 3, 4, 5);
    
    % preallocate seg_idx_l for task-table generation
    seg_idx_l = zeros(64000, 2);
    task_idx = 1;
    t_unit = granularity * vpa(1/sampleRateDAC);
    for idx = (1:length(Y_l) - 1)
        seg_idx_l_added = [];
        Y_spacing = Y_l(idx+1, 1) - Y_l(idx, 2);
        % add y-pulse
        [seg_idx_l(task_idx, 1), seg_idx_l(task_idx, 2)] = deal(y_seg_idx, 1);
        task_idx = task_idx + 1;
        % add x-pulses w/ readout window if spacing exists
        num = fix(Y_spacing/x_seg_len);
        if num > 0
            [seg_idx_l(task_idx, 1), seg_idx_l(task_idx, 2)] = deal(x_seg_idx, num);
            task_idx = task_idx + 1;
        end
        
        % To fill up the remainder - squeeze in x-pulse w/out readout
        % window if there is enough time in between
        rem = mod(Y_spacing, x_seg_len);
        if rem >= x_seg_n_acq_len
            [seg_idx_l(task_idx, 1), seg_idx_l(task_idx, 2)] = deal(x_seg_n_acq_idx, 1);
            task_idx = task_idx + 1;
            rem = rem - x_seg_n_acq_len;
        end
        num_DC_unit = round(rem/t_unit);
        if num_DC_unit > 0
            [seg_idx_l(task_idx, 1), seg_idx_l(task_idx, 2)] = deal(DC_unit_idx, num_DC_unit);
            task_idx = task_idx + 1;
        end
    end
    seg_idx_l = seg_idx_l((1: task_idx - 1), :);
    fprintf("segment index list generated for DTQC.\n");
end