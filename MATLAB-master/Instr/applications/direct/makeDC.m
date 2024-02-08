function [mydcI, mydcQ] = makeDC(length)


segLen = 64*round(length/64); %must be a multiple of 64

max_dac = 2^16-1;
half_dac = floor(max_dac/2);

dacWave = zeros(1, segLen) + half_dac;

mydcI = dacWave;
mydcQ = dacWave;

end