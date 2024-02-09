%% ADC Example
% Genrate and Capture Data with ADC


%% Create WiFi signal

addpath 'C:\Users\Simon\Documents\MATLAB\Examples\R2019b\instrument\BeaconFrameGenerationUsingQCRFSigGenExample'

%SSID = 'TEST_BEACON'; % Network SSID
%beaconInterval = 1; % In Time units (TU)
band = 5;             % Band, 5 or 2.4 GHz
chNum = 140;          % Channel number, corresponds to 5700MHz
Fc = 5.7E+09;

SSID = 'TABOR_AWG'; % Network SSID
beaconInterval = 1; % In Time units (TU)
band = 2.4;         % Band, 5 or 2.4 GHz
chNum = 3;          % Channel number, corresponds to 2422MHz
Fc = 2.422E+09;     % The center frequency of CH 3

% Generate Beacon frame
[mpduBits,fc] = helperGenerateBeaconFrame(chNum, band, beaconInterval, SSID);

cfgNonHT = wlanNonHTConfig;              % Create a wlanNonHTConfig object
cfgNonHT.PSDULength = numel(mpduBits)/8; % Set the PSDU length in bits

% The idle time is the length in seconds of an idle period after each
% generated packet. The idle time is set to the beacon interval.
waveform = wlanWaveformGenerator(mpduBits, cfgNonHT, 'IdleTime', beaconInterval*1024e-6);


%% re-sample to sclk of AWG
Fs = wlanSampleRate(cfgNonHT);           % Get the input sampling rate
sclk = 9e9;
FsNew = sclk/Fs;
waveformReSamp = IqIdealInterpolationWifi (waveform, FsNew);

%% Modulate onto carrier

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

Fc = round(Fc / (sclk / length(carrierWave))) * sclk / length(carrierWave);
% Carrier generation
carrierWave = exp(1i * 2 * pi * Fc * carrierWave);
% Complex carrier multiplied by complex baseband waveform (IQ modulation)
%Modulated signal is just the real part of the complex product
waveformReSamp = real(waveformReSamp .* carrierWave);
waveformReSamp = waveformReSamp.';


%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Users\Simon\Documents\Tabor Electronics\ProteusAwg-10204-Customer\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

%% Create Administrator
% Create instance of the CProteusAdmin class, and open it.
% Note that only a single CProteusAdmin instance can be open, 
% and it must be kept open till the end of the session.

admin = CProteusAdmin(@OnLoggerEvent);

rc = admin.Open();
assert(rc == 0);

%% Access instrument
try
    slotIds = admin.GetSlotIds();
    numSlots = length(slotIds);
    assert(numSlots > 0);
    
    % If there are multiple slots, let the user select one ..
    sId = slotIds(1);
    if numSlots > 1
        fprintf(1, '\n%d slots were found\n', numSlots);
        for n = 1:numSlots
            sId = slotIds(n);
            slotInfo = admin.GetSlotInfo(sId);
            if ~slotInfo.IsSlotInUse
                modelName = slotInfo.ModelName;
                if slotInfo.IsDummySlot
                    fprintf(1, ' * Slot Number: Model %s [Dummy Slot].\n', sId, ModelName);
                else
                    fprintf(1, ' * Slot Number: Model %s.\n', sId, ModelName);
                end
            end
        end
        choice = input('Enter SlotId ');
        fprintf(1, '\n');
        sId = uint32(choice);
    end
    
    % Connect to the selected instrument ..
    should_reset = true;
    inst = admin.OpenInstrument(sId, should_reset);
    instId = inst.InstrId;
    
    % ---------------------------------------------------------------------
    % Send SCPI commands
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    res = inst.SendScpi('*CLS');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':FREQ:RAST 9E9');
    assert(res.ErrCode == 0);    
    
    fprintf('Reset complete\n');

    % ---------------------------------------------------------------------
    % setup waveform memory on CH1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);
    
    sampleRateDAC = sclk;
    segLen = 10240000;
    bits = 8;
    amplitude = 1;
    cycles = 1092;

    % Define segment 1 
    res = inst.SendScpi(':TRAC:DEF 1,10240000');
    assert(res.ErrCode == 0);

    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 
    
    waveformReSampTrunk = waveformReSamp(1:10240000); %truncate so its divisable by 64
    bits=8;
    dacSignal = ampScale(bits, waveformReSampTrunk); 

    res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
    assert(res.ErrCode == 0);
    

    % ---------------------------------------------------------------------
    % Play segment 1 in channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);    

    res = inst.SendScpi(':SOUR:FUNC:SEG 1');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SOUR:VOLT 0.3');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);

    fprintf('Waveform generated and playing\n');

    % ---------------------------------------------------------------------
    % Operate the ADC with direct functions 
    % ---------------------------------------------------------------------


    sampleRateADC = 2e9;
    
    frameLen = 2560032; % integer div of 96
    offLen = 0 
    ch = 1; %CH 1 of Digitizer
    on = 1;
    off = 0;
    numFrames = 1;
    frameNumber = 1;
    
    % ----------
    % ADC Config
    % ----------
    rc = inst.SetAdcCaptureEnable(off);%Disable capture during configuration
    assert(rc == 0);
    
    rc = inst.SetAdcDualChanMode(on);%Dual channel mode for ADC
    assert(rc == 0);
    
    rc = inst.SetAdcSamplingRate(sampleRateADC);%Set sampling rate of ADC
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(0,1000);%Set full scale of channel 1
    assert(rc == 0);

    % The prototype of SetAdcCaptureTrigSource has been changed ..
    trigSource = ch + 1; % self-trigger from ADC channel-1
    trigSource = 0;
    adcChanInd = ch - 1; % ADC Channel 1
    rc = inst.SetAdcCaptureTrigSource(adcChanInd, trigSource);
    assert(rc == 0);
    
    rc = inst.SetAdcSelfTrigThresh(adcChanInd, -0.05);%Set trigger level threshold, CH1, 0.05V
    assert(rc == 0);
    
    rc = inst.SetAdcAcquisitionEn(on,on); %ch1 and ch2 enable
    assert(rc == 0);
    
    rc = inst.SetAdcFramesLayout(numFrames, frameLen);
    assert(rc == 0);
    
    % ---------------------------------------------------------------------
    % Capture Parmeters
    % ---------------------------------------------------------------------
       
    %rc = inst.SetAdcMultiFrameModeEn(on);
    %assert(rc == 0);
     
    % ---------------------------------------------------------------------
    % ADC Capture
    % ---------------------------------------------------------------------
              


    % allocate NET array
    netArray = NET.createArray('System.UInt16', frameLen * numFrames);

    fprintf('Capture set\n');
    h = figure(1);
    forever = 1; 
    while forever 
        drawnow
        isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
        if isKeyPressed
            break
        end

        rc = inst.SetAdcCaptureEnable(on);
        assert(rc == 0);
        
        % Generate software-trigger for ADC capturing
        rc = inst.GenerateAdcSoftTrig();
        assert(rc == 0);        

        % Wait till the capture completes
        status = inst.ReadAdcCaptureStatus();
        for i = 1 : 2500
            if status ~= 0
                break;
            end
            status = inst.ReadAdcCaptureStatus();
        end
        
        rc = inst.SetAdcCaptureEnable(off);
        assert(rc == 0);
        
        % read the captured frame
        rc = inst.ReadMultipleAdcFrames(adcChanInd, frameNumber, numFrames, netArray);

        samples = int16(netArray);
        
        maxSig = max(samples);
        minSig = min(samples);
        
        signalRangePkPk = (minSig*-1)+maxSig;

        meanSig = mean(samples);
        dataReadTimeDC=samples-meanSig;
       
        n = 1;
        while n < frameLen
            
            if dataReadTimeDC(n) > 200
                dataReadTimeDC(n);
                break;
            end
            n = n+1;
        end
        burstLen = 240000;
        dataReadTimeDCGated = dataReadTimeDC(n:n+burstLen);
        
        X = fftshift(fft(dataReadTimeDCGated));
        N=burstLen+1;
        dF = sampleRateADC/N;                      % hertz
        f = -sampleRateADC/2:dF:sampleRateADC/2-dF;
        
        startPointFull = int32(N/2);
        stopPointFull = int32(length(f));
        
        underSampFreq = Fc-sampleRateADC; % ~400MHz
        span = Fs*2; % 40MHz
        hzPerPoint = sampleRateADC / N;
        startFreq = underSampFreq - (Fs/2);
        stopFreq = underSampFreq + (Fs/2);
        
        startPointsZoom = startFreq/hzPerPoint;
        stopPointsZoom = stopFreq/hzPerPoint;
        
        tiledlayout(2,2);
        nexttile;
        plot(samples);
        title('Captured Waveform');
        xlabel('Points')
        nexttile;
        plot(dataReadTimeDCGated);
        title('Gated Waveform');
        xlabel('Points');
        nexttile;
        plot(f(startPointFull:stopPointFull),abs(X(startPointFull:stopPointFull))/N);
        title('Spectrum Full Span 0Hz-1GHz');
        xlabel('Hz');
        nexttile;
        plot(f(startPointFull+startPointsZoom:startPointFull+stopPointsZoom),abs(X(startPointFull+startPointsZoom:startPointFull+stopPointsZoom))/N);
        title('Spectrum Zoom');
        xlabel('Hz');

        %pause(0.05);
        %forever = 0;
    end 
    
    	% deallocate the NET array
    delete(netArray);
    
    %figure(1); clf;
    %samples = samples - 2048;
    %samples = myNormalization(samples, 1);

    % Free the memory space that was allocated for ADC capture
    rc = inst.FreeAdcReservedSpace();
    assert(rc == 0);


% ---------------------------------------------------------------------
% End of the example
% ---------------------------------------------------------------------
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd of Example - %s\n', netStrToStr(res.RespStr));
    close all; % Close all figures
    % It is recommended to disconnect from instrumet at the end
    rc = admin.CloseInstrument(instId);    
    
    % Close the administrator at the end ..
    admin.Close();
catch ME
    admin.Close();
    rethrow(ME)
end

% Function netStrToStr
function str = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch        
        str = '';
    end
end

%% Function - On Logger-Event
function OnLoggerEvent(~, e)
    try       
        % print only:
        % TaborElec.Proteus.CLI.LogLevel.Critical (1)
        % TaborElec.Proteus.CLI.LogLevel.Error    (2)
        % TaborElec.Proteus.CLI.LogLevel.Warning  (3)
        if int32(e.Level) >= 1 &&  int32(e.Level) <= 3
            msg = netStrToStr(e.Message);
            if strlength(msg) > 0
                fprintf(2, '\n ~ %s\n', msg(1));
            end
        end
        System.Diagnostics.Trace.WriteLine(e.Message);
    catch ME
       rethrow(ME) 
    end
end

%% Function - Interpolate
function retval = IqIdealInterpolationWifi (myArray, xFactor)
  
  %expansion by zero-padding
  retval = zeros(1, xFactor * length(myArray));
  retval([1:xFactor:end]) = myArray;
  % "Ideal" Interpolation filter
  lenSinc = 40;  
  mySinc = -lenSinc * xFactor : 1 : lenSinc * xFactor ;  
  mySinc = sinc(mySinc / xFactor);
  myWindow = blackman(length(mySinc));
  myWindow = myWindow.'; 
  mySinc = mySinc .* myWindow;
  %convolution
  retval = cconv(retval, mySinc, length(retval));
  %retval = real(retval);
  
end

%% Function - scale
function dacSignal = ampScale(bits, rawSignal)
 
  maxSig = max(rawSignal);

  verticalScale = ((2^bits)/2)-1;

  vertScaled = (rawSignal / maxSig) * verticalScale;
  dacSignal = uint8(vertScaled + verticalScale);
  dacSignal = typecast(dacSignal, 'uint8');

  if bits > 8
      dacSignal16 = [];
      sigLen = length(dacSignal);
      k=1;
      for j = 1:2:sigLen*2;
        dacSignal16(j) = bitand(dacSignal(k), 255);
        dacSignal16(j+1) = bitshift(dacSignal(k),-8);
        k = k + 1;
      end
      dacSignal16;
      dacSignal = dacSignal16;
  end
  
end



