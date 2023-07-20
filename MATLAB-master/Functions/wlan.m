% Generated by MATLAB(R) 9.7 (R2019b) and WLAN Toolbox 2.2 (R2019b).
% Generated on: 26-Feb-2020 15:21:31

%% Generate 802.11ax Waveform
% 802.11ax configuration:
heSUCfg = wlanHESUConfig('ChannelBandwidth', 'CBW20', ...
    'NumTransmitAntennas', 1, ...
    'NumSpaceTimeStreams', 1, ...
    'SpatialMapping', 'Direct', ...
    'PreHESpatialMapping', false, ...
    'MCS', 4, ...
    'DCM', false, ...
    'ChannelCoding', 'LDPC', ...
    'APEPLength', 100, ...
    'GuardInterval', 3.2, ...
    'HELTFType', 4, ...
    'UplinkIndication', false, ...
    'BSSColor', 0, ...
    'SpatialReuse', 0, ...
    'TXOPDuration', 127, ...
    'HighDoppler', false);


numPackets = 1;
% input bit source:
in = randi([0, 1], 1000, 1);


% waveform generation:
waveform = wlanWaveformGenerator(in, heSUCfg, ...
    'NumPackets', numPackets, ...
    'IdleTime', 0, ...
    'ScramblerInitialization', 93, ...
    'WindowTransitionTime', 1e-07);

Fs = wlanSampleRate(heSUCfg); 								 % sample rate of waveform

Fc = 0.9E+09;
sclk = 2.5e9;
FsNew = sclk/Fs;

waveformReSamp = IqIdealInterpolation (waveform(1:1664), FsNew);
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

%% Visualize 802.11ax Waveform
% Spectrum Analyzer
spectrum = dsp.SpectrumAnalyzer('SampleRate', Fs);
spectrum(waveform);
release(spectrum);

% RU Assignment and Allocated Subcarriers:
showAllocation(heSUCfg);

