SSID = 'TEST_BEACON'; % Network SSID
beaconInterval = 1; % In Time units (TU)
band = 5;             % Band, 5 or 2.4 GHz
chNum = 140;          % Channel number, corresponds to 5700MHz

% Generate Beacon frame
[mpduBits,fc] = helperGenerateBeaconFrame(chNum, band, beaconInterval, SSID);

cfgNonHT = wlanNonHTConfig;              % Create a wlanNonHTConfig object
cfgNonHT.PSDULength = numel(mpduBits)/8; % Set the PSDU length in bits

% The idle time is the length in seconds of an idle period after each
% generated packet. The idle time is set to the beacon interval.
waveform = wlanWaveformGenerator(mpduBits, cfgNonHT, 'IdleTime', beaconInterval*1024e-6);
Fs = wlanSampleRate(cfgNonHT);           % Get the input sampling rate