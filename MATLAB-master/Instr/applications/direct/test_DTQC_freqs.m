% Some numbers that work
sampleRateDAC = 675000000;
granularity = 64;
gamma_len = 100e-6;

% Assuming tacq = 12e-6
gamma_pulse_len = gamma_len + 36e-6;
T1_l = (2e-3: 1e-5: 7e-3);
T2_l = (2e-3: 1e-5: 6e-3);

%initialize gamma1 array and gamma2 array

for idx1 = (1:length(T1_l))
    for idx2 = (1:length(T2_l))
        gamma1_l = [0,gamma_pulse_len];
        gamma2_l = [gamma_pulse_len,2*gamma_pulse_len];
        T1 = T1_l(idx1);
        T2 = T2_l(idx2);
        [T1, T2] = deal(T1, T2);
        time_arr = round_to_DAC_freq([T1, T2], sampleRateDAC, granularity);
        [T1, T2] = deal(time_arr(1), time_arr(2));
        [freq1, freq2] = deal(1/T1, 1/T2);
        seq_time = 5;
        
        % add all possible gamma1 pulses
        while gamma1_l(end, 2) + T1 < seq_time
            gamma1_l = [gamma1_l; gamma1_l(end, 1)+T1, gamma1_l(end, 2)+T1];
        end

        % add all possible gamma2 pulses
        while gamma2_l(end, 2) + T2 < seq_time
            gamma2_l = [gamma2_l; gamma2_l(end, 1)+T2, gamma2_l(end, 2)+T2];
        end
        gamma_l = [gamma1_l; gamma2_l];
        gamma_l = sortrows(gamma_l, 1);
        % fprintf("length of gamma_l: %d \n", length(gamma_l));
        % need to check if the any two gamma pulses overlap
        overlap = false;
        for idx = (1: length(gamma_l)-1)
            if gamma_l(idx, 2) > gamma_l(idx+1, 1)
                overlap = true;
                % fprintf("%d, %d \n", gamma_l(idx, 2), gamma_l(idx+1, 1));
            end
        end
        if overlap == false
            fprintf("overlap is false\n")
            fprintf('%d, %d \n', T1, T2);
        end
    end
end

fprintf('%s \n', mat2str(overlap));
% LOGIC
% if I can fit a spin-lock without readout window -- do it
% if length is not enough, then don't but extend the last pulse.
% every last pulse will have a readout window.

fprintf("Done \n");