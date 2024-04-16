function rounded_t_arr = round_to_DAC_freq(time_arr, sampleRateDAC, granularity)
    rounded_t_arr = [];
    for idx = (1: length(time_arr))
        t_unit = granularity * vpa(1/sampleRateDAC);
        t = round(time_arr(idx)/t_unit)*t_unit;
        rounded_t_arr(end+1) = t;
    end
end