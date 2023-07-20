%% 802.11 OFDM Beacon Frame Generation and Transmission Using Quick-Control RF Signal Generator
% This example shows how to generate packets containing MAC beacon frames
% suitable for baseband simulation or over-the-air transmission using WLAN
% Toolbox&trade;, Instrument Control Toolbox&trade; and Keysight Technologies&reg; RF signal generator.

% Copyright 2017-2018 The MathWorks, Inc. 

%% Introduction
% In this example WLAN Toolbox is used to create an IEEE&reg;
% 802.11&trade; beacon frame. Using Instrument Control Toolbox, the generated beacon frame is downloaded to Keysight 
% Technologies N517B signal generator for over-the-air transmission.
% Beacon frame is a type of management frame that identifies a basic 
% service set (BSS) formed by a number of 802.11 devices. The access point 
% of a BSS periodically transmits the beacon frame to establish and maintain the network. 
% A WiFi device can be used to view this beacon frame transmitted by the RF 
% Signal Generator.

%%
% <<../BeaconFrameGenerationRx.png>>
%
% For more information on beacon 
% frame generation using WLAN Toolbox, please refer 
% <https://www.mathworks.com/help/wlan/examples/802-11-ofdm-beacon-frame-generation.html
% 802.11 OFDM Beacon Frame Generation.>

%% Requirements
% To run this example you need:
% 
% * Keysight Technologies N5172B signal generator
% * Keysight VISA version 17.3
% * IVI-C driver for Keysight Technologies N5172B signal generator
% * National Instruments&trade; IVI&reg; compliance package version 16.0.1.2 or higher
% * WLAN Toolbox
% * Instrument Control Toolbox

%% Create IEEE 802.11 Beacon Frame
% The beacon packets are periodically transmitted with a beacon interval of 100 TU, where
% 1 TU represents 1024 microseconds time interval between successive
% beacons.The MAC frame bits for the beacon frames are generated using the helper
% function <matlab:edit('helperGenerateBeaconFrame.m')
% helperGenerateBeaconFrame>.

SSID = 'TEST_BEACON'; % Network SSID
beaconInterval = 100; % In Time units (TU)
band = 5;             % Band, 5 or 2.4 GHz
chNum = 52;           % Channel number, corresponds to 5260MHz

% Generate Beacon frame
[mpduBits,fc] = helperGenerateBeaconFrame(chNum, band, beaconInterval, SSID);

%% Create IEEE 802.11 Beacon Packet
% A beacon packet is synthesized using <matlab:doc('wlanWaveformGenerator')
% wlanWaveformGenerator> with a non-HT format configuration object. In this
% example an object is configured to generate a beacon packet of 20 MHz
% bandwidth, 1 transmit antenna and BPSK rate 1/2 (MCS 1).

cfgNonHT = wlanNonHTConfig;              % Create a wlanNonHTConfig object
cfgNonHT.PSDULength = numel(mpduBits)/8; % Set the PSDU length in bits

% The idle time is the length in seconds of an idle period after each
% generated packet. The idle time is set to the beacon interval.
txWaveform = wlanWaveformGenerator(mpduBits, cfgNonHT, 'IdleTime', beaconInterval*1024e-6);
Rs = wlanSampleRate(cfgNonHT);           % Get the input sampling rate

%% Create a RF Signal Generator Object
% Quick-Control RF Signal Generator is used to download and transmit the
% baseband waveform, |txWaveform|, generated by WLAN Toolbox.
rf = rfsiggen();
%%
% Discover all the available instrument resources you can connect to, using the |resources| method.
rf.resources
%%
% Discover all the available instrument drivers, using |drivers| method.
rf.drivers

%% Connect to Signal Generator
% Set |Resource| and |Driver| properties before connecting to the object. The IP address of Keysight Technologies N5172B signal generator 
% is _172.28.21.217_ , hence the resource specified will be 'TCPIP0::172.28.21.217::inst0::INSTR'
rf.Resource = 'TCPIP0::172.28.21.217::inst0::INSTR';
rf.Driver = 'AgRfSigGen';
% Connect to the instrument
connect(rf);

%% Download Waveform
% Download the waveform, |txWaveform|, to the instrument with sampling rate
% |Rs|.
download(rf, transpose(txWaveform), Rs);

%% Transmit the Waveform
% Call |start| to start transmitting waveform using specified 
% centerFrequency, outputPower and loopCount. 
centerFrequency = fc;
outputPower = 0;
loopCount = Inf;
start(rf, centerFrequency, outputPower, loopCount);

%%
% Once the signal generator is transmitting the beacon, you can test by
% scanning for wireless network using a Wi-Fi device. You should now see a
% _TEST_BEACON_ SSID in the list of available networks.

%% Clean Up
% When you have finished transmitting, stop the waveform output, disconnect the 
%  |rfsiggen| object from the signal generator, and remove it from the workspace.
stop(rf);
disconnect(rf);
clear rf


%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperGenerateBeaconFrame.m') helperGenerateBeaconFrame.m>
% * <matlab:edit('helperWLANChannelFrequency.m') helperWLANChannelFrequency.m>
  
