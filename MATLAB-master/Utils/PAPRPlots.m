magWaveform = abs(waveform);

%Time Gate function
gatePoint = length(magWaveform); 
packetOff = zeros(100,1);
packetOff = packetOff.';
n=1;
while n < length(magWaveform)
    if magWaveform(n:n+99)== packetOff
        gatePoint = n;
        break
    end
    n = n + 1;
end

gatedWaveform =  magWaveform(1:gatePoint);

% Calculate CCDF
[xCCDF, yCCDF] = CCDF(gatedWaveform, 100);

% Calculate Peak to Avg
envelopeWfm = gatedWaveform .* gatedWaveform;
papr = max(envelopeWfm) / mean(envelopeWfm);
papr = 10 * log10(papr);


%PLOTS
tiledlayout(2,2);
nexttile;
plot(real(waveform));
title('Waveform');
nexttile;
plot(abs(fft(waveformReSamp)))
title('FFT Waveform');
nexttile;
plot(gatedWaveform);
title('Gated Magnitude');
nexttile;
semilogy (xCCDF, yCCDF);
title(strcat(waveFormDescription, ' CCDF Plot, PAPR = ', num2str(papr), 'dB'));


spectrum = dsp.SpectrumAnalyzer('SampleRate', sclk);
spectrum(waveformReSamp);
release(spectrum);