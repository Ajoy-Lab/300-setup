SSID = 'PR'; % Network SSID
beaconInterval = 1; % In Time units (TU)
band = 5;             % Band, 5 or 2.4 GHz
chNum = 140;          % Channel number, corresponds to 5700MHz
Fc = 1E+09;

% Generate Beacon frame
[mpduBits,fc] = helperGenerateBeaconFrame(chNum, band, beaconInterval, SSID);

cfgNonHT = wlanNonHTConfig;              % Create a wlanNonHTConfig object
cfgNonHT.PSDULength = numel(mpduBits)/8; % Set the PSDU length in bits

% The idle time is the length in seconds of an idle period after each
% generated packet. The idle time is set to the beacon interval.
waveform = wlanWaveformGenerator(mpduBits, cfgNonHT, 'IdleTime', beaconInterval*1024e-6);
Fs = wlanSampleRate(cfgNonHT);           % Get the input sampling rate
sclk = 2.5e9;
FsNew = sclk/Fs;

waveformReSamp = IqIdealInterpolation (waveform(1:22848), FsNew);
% Carrier Waveform creation
carrierWave = 0:(length(waveformReSamp) - 1);
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
% Carrier generation
carrierWave = exp(1i * 2 * pi * Fc * carrierWave);
% Complex carrier multiplied by complex baseband waveform (IQ modulation)
%Modulated signal is just the real part of the complex product
waveformReSamp = real(waveformReSamp .* carrierWave);
%plot(abs(fft(waveformReSamp(1:100000))))
waveformReSamp = waveformReSamp.';