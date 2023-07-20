function [mpduBits,fc] = helperGenerateBeaconFrame(chNumber,band,beaconInterval,ssid)
%helperGenerateBeaconFrame Generate MAC beacon frame
%
%   [MPDUBITS,FC] = helperGenerateBeaconFrame(CHNUMBER,BEACONINTERVAL, ...
%   SSID) returns the beacon frame bits and carrier frequency for the
%   specified chNumber, beaconInterval and SSID.

% Copyright 2016 The MathWorks, Inc.

%#codegen

% Get carrier frequency for the specified channel number and band
fc = helperWLANChannelFrequency(chNumber, band);
            
supportedRates = [6 9 12 18 24 36 48 54];
sendCFParameters = false;
numMPDUOctets = 24 ...                                 % MPDU Header
               + 8 ...                                 % Timestamp
               + 2 ...                                 % Beacon interval
               + 2 ...                                 % Capability
               + 2 + length(uint8(ssid)) ...           % SSID Element
               + 2 + length(supportedRates) ...        % Supported supportedRates Element
               + 3 + ...                               % DS Parameter set
               + 6 ;                                   % TIM parameters

baseSeqNum = 1395; % Fixed to 1395
            
% Generate FCS and from MPDU, IEEE Std 802.11-2012 Section 8.2.4
fcsGen = comm.CRCGenerator([32 26 23 22 16 12 11 10 8 7 5 4 2 1 0], ...
            'InitialConditions', 1, ...
            'DirectMethod', true, ...
            'FinalXOR', 1);

% Set output
p80211 = struct( ...
        'BeaconInterval',   beaconInterval, ...
        'SSID',             uint8(ssid), ...
        'SupportedRates',   supportedRates, ...
        'ChannelNumber',    chNumber, ...
        'SendCFParameters', sendCFParameters, ...
        'NumMPDUOctets',    numMPDUOctets, ...
        'TSF',              zeros(2, 1, 'uint32'), ...
        'TSFD',             beaconInterval*1024, ...
        'SequenceNumber',   baseSeqNum, ...
        'FrameType',        'Beacon');
            
% Get MAC frame
mpduBits = generateMPDU(p80211, fcsGen);
            
end
       
function y = generateMPDU(p80211,fcsGen)
        
    p80211.TSF = p80211.TSF + p80211.TSFD; % Increment TSF
    p80211.SequenceNumber = p80211.SequenceNumber + 1; % Increment sequence number

    [mpduHeader, mpduBody] = getNonHTOFDMBeacon(p80211);
    payload = serializeMPDUData(mpduHeader, mpduBody);

    % Data in payload is octets, convert to bits ordered by octet (right bit MSB)
    bits = de2bi(payload).';
    y = step(fcsGen,double(bits(:)));  % Append FCS bits
          
end

function [H, B] = getNonHTOFDMBeacon(p80211)
%   [H,B] = getNonHTOFDMBeacon(p80211) returns a beacon frame header and
%   body fields for the 802.11 WLAN beacon. H and B are structures that
%   contain the beacon frame header and body fields. Input p80211
%   represents the parameter structure containting beacon information
%   fields within a beacon frame body as defined in IEEE Std 802.11-2007
%   Section 7.2.3.1.

type    = 0;    % Management frame (00)
subtype = 8;    % Beacon frame (1000)

% MPDU header
H = getMPDUHeader(type, subtype, p80211.SequenceNumber);

% Frame body
B = getBeaconMPDUBody(p80211.TSF, p80211);
end

%--------------------------------------------------------------------------
function mpduHeader = getMPDUHeader(type, subtype, seqNum)
    mpduHeader.FrameCtrl       = getFrameControl(type, subtype);
    mpduHeader.DurationID      = getDurationID;
    mpduHeader.Address1        = [255 255 255 255 255 255]; % Dummy destination Addrs
    mpduHeader.Address2        = [200 96 0 147 241 221];    % Dummy station Addrs
    mpduHeader.Address3        = [200 96 0 147 241 221];    % Dummy BSSID
    mpduHeader.SequenceControl = getSequenceControl(seqNum);
    mpduHeader.Address4        = [0 0 0 0 0 0];             % Dummy value
end

% Get Frame Control
function frameCtrl = getFrameControl(type, subtype)
    frameCtrl.ProtocolVersion = uint8(0);
    frameCtrl.Type            = uint8(type);
    frameCtrl.Subtype         = uint8(subtype);
    frameCtrl.ToDS            = uint8(0);
    frameCtrl.FromDS          = uint8(0);
    frameCtrl.MoreFragments   = uint8(0);
    frameCtrl.Retry           = uint8(0);
    frameCtrl.PowerManagement = uint8(0);
    frameCtrl.MoreData        = uint8(0);
    frameCtrl.ProtectedFrame  = uint8(0);
    frameCtrl.Order           = uint8(0);
end

% Get Duration
function durationIDOct = getDurationID()
    durationID = [de2bi(2365,15) 0]; % Duration value
    durationIDOct = convertToOctets(durationID);
end

% Get Sequence number
function sequenceControl = getSequenceControl(seqNum)
    fragNum = uint16(0);
    seqNum = uint16(seqNum);
    sequenceControl = convertToOctets(de2bi(bitor(bitshift(seqNum,4), fragNum), 16)); 
end

% Get beacon body fields
function beaconMPDUBody = getBeaconMPDUBody(TSF, p80211)
    BITS_PER_OCTET = 8;

    % 1 - Timestamp
    beaconMPDUBody.TimeStamp = zeros(1, 8, 'uint8');
    for p=0:3
      beaconMPDUBody.TimeStamp(p+1) = ...
        bitand(bitshift(TSF(1), -p*BITS_PER_OCTET),255);
    end
    for p=0:3
      beaconMPDUBody.TimeStamp(p+5) = ...
        bitand(bitshift(TSF(2), -p*BITS_PER_OCTET),255);
    end
      
    % 2 - Beacon interval
    beaconMPDUBody.BeaconInterval = dec2int8(uint16(p80211.BeaconInterval),2);

    % 3 - Capability
    beaconMPDUBody.Capability = getCapabilities();

    % 4 - Info elements
    infoElementStruct = struct('ID', uint8(0), ...
      'Length', uint8(0), ...
      'Value', zeros(1,32,'uint8'));
    infoElemCount = 1;

    % 5 - SSID
    beaconMPDUBody.InfoElements = repmat(infoElementStruct, 1, 25);
    beaconMPDUBody.InfoElements(infoElemCount) = ...
      getSSIDElement(infoElementStruct, p80211.SSID);
    infoElemCount = infoElemCount + 1;

    % 6 - Supported rates
    beaconMPDUBody.InfoElements(infoElemCount) = ...
      getSupportedDataRates(infoElementStruct, p80211.SupportedRates);
    infoElemCount = infoElemCount + 1;

    % 7 - Channel number
    beaconMPDUBody.InfoElements(infoElemCount) = ...
      getDSParameterSet(infoElementStruct, p80211.ChannelNumber);
    infoElemCount = infoElemCount + 1;

    % 8 - Traffic indication MAP
    beaconMPDUBody.InfoElements(infoElemCount) = ...
      getTIMParameterSet(infoElementStruct);

    beaconMPDUBody.NumInfoElements = infoElemCount;
end

function capability = getCapabilities()
    capability.ESS                = true;
    capability.IBSS               = false;
    capability.CFPollable         = false;
    capability.CFPollRequest      = false;
    capability.Privacy            = false;
    capability.ShortPreamble      = false;
    capability.PBCC               = false;
    capability.ChannelAgility     = false;
    capability.SpectrumManagement = false;
    capability.QoS                = false;
    capability.ShortSlotTime      = false;
    capability.APSD               = false;
    capability.Reserved           = false;
    capability.DSSOFDM            = false;
    capability.DelayedBlockAck    = false;
    capability.ImmediateBlockAck  = false;
end

% SSID
function element = getSSIDElement(element, ssid)
    len = length(ssid);
    element.ID = uint8(0);
    element.Length = uint8(len);
    element.Value(1:len) = ssid;
end

% Supported rates
function element = getSupportedDataRates(element, supportedRates)
    len = length(supportedRates);
    element.ID = uint8(1);
    element.Length = uint8(len);
    % Increments of 0.5 Mbps with MSB set to 1
    element.Value(1:len) = uint8(ceil(supportedRates / 0.5) + 2^7);
end

% Channel
function element = getDSParameterSet(element, chanNum)
    element.ID = uint8(3);
    element.Length = uint8(1);
    element.Value(1) = uint8(chanNum);
end

% TIM
function element = getTIMParameterSet(element)
    element.ID = uint8(5);
    element.Length = uint8(4);
    element.Value(1) = uint8(0);
    element.Value(2) = uint8(0);
    element.Value(3) = uint8(0);
    element.Value(4) = uint8(0);
end

function out = convertToOctets(in)
    BITS_PER_OCTET = 8;

    numOctets = length(in)/BITS_PER_OCTET;
    mask = 2.^(0:BITS_PER_OCTET-1)';

    out = zeros(1,numOctets,'uint8');
    for p=1:numOctets
      out(p) = double(in((p-1)*BITS_PER_OCTET+1:p*BITS_PER_OCTET))*mask;
    end
end

function out = dec2int8(in, numInt8)
    BITS_PER_OCTET = 8;

    out = zeros(1, numInt8, 'uint8');
    for p=1:numInt8
      out(p) = bitand(in, 255);
      in = bitshift(in, -BITS_PER_OCTET);
    end
end

function [D, L] = serializeMPDUData(mpduHeader, mpduBody)
%SERIALIZEMPDUDATA Serialize the MPDU data in a uint8 array
%   [D,L] = serializeMPDUData(H,B) serializes the MPDU data stored in
%   the header structure H and body structure B and then returns the data
%   D, which is a uint8 array. L represents the number of bits in D.

MaximumMPDULength = 68; % length in octets of beacon information fields
D = zeros(MaximumMPDULength, 1, 'uint8'); % + 32-bit CRC

% Process MPDU header
p=uint16(1);
D(p) = mpduHeader.FrameCtrl.ProtocolVersion ...
    + bitshift(mpduHeader.FrameCtrl.Type, 2) ...
    + bitshift(mpduHeader.FrameCtrl.Subtype, 4); p=p+1;
D(p) = mpduHeader.FrameCtrl.ToDS ...
    + bitshift(mpduHeader.FrameCtrl.FromDS, 1) ...
    + bitshift(mpduHeader.FrameCtrl.MoreFragments, 2) ...
    + bitshift(mpduHeader.FrameCtrl.Retry, 3) ...
    + bitshift(mpduHeader.FrameCtrl.PowerManagement, 4) ...
    + bitshift(mpduHeader.FrameCtrl.MoreData, 5) ...
    + bitshift(mpduHeader.FrameCtrl.ProtectedFrame, 6) ...
    + bitshift(mpduHeader.FrameCtrl.Order, 7); p=p+1; 

D(p:p+1) = mpduHeader.DurationID; p=p+2;
D(p:p+5) = mpduHeader.Address1; p=p+6;
D(p:p+5) = mpduHeader.Address2; p=p+6;
D(p:p+5) = mpduHeader.Address3; p=p+6;
D(p:p+1) = mpduHeader.SequenceControl; p=p+2;

% Process MPDU payload (Beacon frame)
D(p:p+7) = mpduBody.TimeStamp; p=p+8;
D(p:p+1) = mpduBody.BeaconInterval; p=p+2;
D(p) = uint8(mpduBody.Capability.ESS) ...
    + bitshift(uint8(mpduBody.Capability.IBSS), 1) ...
    + bitshift(uint8(mpduBody.Capability.CFPollable), 2) ...
    + bitshift(uint8(mpduBody.Capability.CFPollRequest), 3) ...
    + bitshift(uint8(mpduBody.Capability.Privacy), 4) ...
    + bitshift(uint8(mpduBody.Capability.ShortPreamble), 5) ...
    + bitshift(uint8(mpduBody.Capability.PBCC), 6) ...
    + bitshift(uint8(mpduBody.Capability.ChannelAgility), 7); p=p+1;
D(p) = uint8(mpduBody.Capability.SpectrumManagement) ...
    + bitshift(uint8(mpduBody.Capability.QoS), 1) ...
    + bitshift(uint8(mpduBody.Capability.ShortSlotTime), 2) ...
    + bitshift(uint8(mpduBody.Capability.APSD), 3) ...
    + bitshift(uint8(mpduBody.Capability.Reserved), 4) ...
    + bitshift(uint8(mpduBody.Capability.DSSOFDM), 5) ...
    + bitshift(uint8(mpduBody.Capability.DelayedBlockAck), 6) ...
    + bitshift(uint8(mpduBody.Capability.ImmediateBlockAck), 7); p=p+1;

% Add Information elements
for elemCount = 1:mpduBody.NumInfoElements
    D(p) = mpduBody.InfoElements(elemCount).ID; p=p+1;
    len = uint16(mpduBody.InfoElements(elemCount).Length);
    D(p) = mpduBody.InfoElements(elemCount).Length; p=p+1;
    D(p:p+len-1) = ...
      mpduBody.InfoElements(elemCount).Value(1:len); p=p+len;
end

L = p - 1;

BITS_PER_OCTET = 8;
L = L * BITS_PER_OCTET;

end

% [EOF]