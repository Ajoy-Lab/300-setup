
FsNew = round(sclk/Fs);

envelopeWfm = abs(waveform);
envelopeWfm = envelopeWfm .* envelopeWfm;

waveformInterp = IqIdealInterpolation (waveform, FsNew);
% Carrier Waveform creation
carrierWave = 0:(length(waveformInterp) - 1);
carrierWave = carrierWave ./ sclk;

% second Nyquist band generation
if Fc > sclk / 2
    Fc = sclk - Fc;
    % one way to reverse the spectrum is changing the sign of the time so
    % carrier rotation in the complex plane goes in the opposite direction
    carrierWave = -carrierWave;
end

%carrierWave = complex(carrierWave);
% Integer number of cycles
Fc = round(Fc / (sclk / length(carrierWave))) * sclk / length(carrierWave);

waveformUpCon = IqModulator2(waveformInterp, Fc, 0, sclk);
% Complex carrier multiplied by complex baseband waveform (IQ modulation)
%Modulated signal is just the real part of the complex product
%waveformReSamp = real(waveformReSamp .* carrierWave);
%plot(abs(fft(waveformReSamp(1:100000))))
waveformReSamp = waveformUpCon.';

