

%% ADC Example
% Capture Data with ADC

clear;
close;


%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

admin = CProteusAdmin(@OnLoggerEvent);
rc = admin.Open();
assert(rc == 0);


%% Use the administrator, and close it at the end
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
    %% 
    
    % ---------------------------------------------------------------------
    % Send SCPI commands
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    res = inst.SendScpi('*CLS'); % clear
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST'); % reset
    assert(res.ErrCode == 0);
    
    sampleRateDAC =9e9;
    res = inst.SendScpi(':FREQ:RAST 9E9'); % set sample clock
    assert(res.ErrCode == 0);    
    
    fprintf('Reset complete\n');

    % ---------------------------------------------------------------------
    % Download sine waveform of 1024 points to segment 1 of channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1'); % select channel 1
    assert(res.ErrCode == 0);
    
    bits = 8;
    
    rampTime = 1/750;
    fStart = 363e6;
    fStop = 367e6;
    dt = 1/sampleRateDAC;
    dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
    fprintf('waveform length - ');
    fprintf(num2str(length(dacSignal)));
    fprintf('\n') ;
    
    % Define segment 1 
    res = inst.SendScpi(strcat(':TRAC:DEF 1,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
    assert(res.ErrCode == 0);

    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 

    % Download the binary data to segment 1
    res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SOUR:CFR 3.4E9');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SOUR:MODE IQM');
    assert(res.ErrCode == 0);
    
    segLen = 40960;
    segLenDC = 50048;
    cycles = 800;
    
    dacSignalMod = sine(cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)
    
    dacSignalDC = [];
    for i = 1 : segLenDC
      dacSignalDC(i) = 127;
    end
    
    dacSignal = [dacSignalMod dacSignalDC];

     
%     %Set Generator to Trigger mode
%     res = inst.SendScpi(':INIT:CONT OFF');
%     assert(res.ErrCode == 0);
% 
%     %Set trigger source to TRG1
%     res = inst.SendScpi(':TRIG:SOUR:ENAB TRG1');
%     assert(res.ErrCode == 0);
% 
%     %Select trigger to program
%     res = inst.SendScpi(':TRIG:SEL TRG1');
%     assert(res.ErrCode == 0);
% 
%     %Enable TRG1
%     res = inst.SendScpi(':TRIG:STAT ON');
%     assert(res.ErrCode == 0);

%     %Set trigger level
%     res = inst.SendScpi(':TRIG:LEV 0.2');
%     assert(res.ErrCode == 0);

    
    % ---------------------------------------------------------------------
    % Play segment 1 in channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);    

    res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 1');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':SOUR:VOLT 0.6');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);

    fprintf('Waveform generated and playing\n');

    % ---------------------------------------------------------------------
    % Operate the ADC with direct functions 
    % ---------------------------------------------------------------------

    sampleRateADC = 1e9;

    readLen = 33888; % integer div of 96
    readLenSeconds = readLen/1e9;
    chanIndex = 0; %CH 1 of Digitizer
    off = 0;
    on = 1;
    numberOfFrames = 1;
    frameNumber = 1;
    netArray = NET.createArray('System.UInt16', readLen);
    
    %netArray=zeros(1,readLen);

    % Turn on ADC dual-channels mode (state = 1)
    rc = inst.SetAdcDualChanMode(on);
    assert(rc == 0);

    rc = inst.SetAdcSamplingRate(sampleRateADC);
    assert(rc == 0);
    
    rc = inst.SetAdcDualChanMode(off);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(chanIndex,1000); % mv - 500, 800, 1000
    assert(rc == 0);
    
    rc = inst.SetAdcFramesLayout(numberOfFrames, readLen);
    assert(rc == 0);
    
    fprintf('Capture size set\n');
 
    fprintf('Trigger set - ');
    trigSource = 2; % 2 = Self, 1 = External
    rc = inst.SetAdcCaptureTrigSource(chanIndex, trigSource); % 1 means external trigger, 2 means self trigger
    assert(rc == 0); 

    if trigSource == 1
        rc = inst.SetAdcExternTrigPattern(0);
        assert(rc==0) 
        rc = inst.SetAdcExternTrigGateModeEn(on);
        assert(rc==0);
        rc = inst.SetAdcExternTrigThresh(0,0.5);
        assert(rc==0);
        rc = inst.SetAdcExternTrigDelay(0.0000025);
        assert(rc==0); 
        fprintf('External\n');
    end

    if trigSource == 2
        %rc = inst.SetAdcSelfTrigThresh(0,0.5); % does not work
        %assert(rc==0);
        fprintf('Self\n');
    end
   
    rc = inst.SetAdcAcquisitionEn(on,off);
    assert(rc == 0);
   
    fprintf('Aquisition Enabled, starting meas loop......\n');
    fprintf('Hit space to exit\n');
    
    h = figure(1);
    forever = 1; 
    while forever 
        drawnow
        isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
        if isKeyPressed
            break
        end
        
        % Wait until the capture completes
        status = inst.ReadAdcCaptureStatus();
        for i = 1 : 250
            if status ~= 0
                break;
            end
            pause(0.01);
            status = inst.ReadAdcCaptureStatus();
        end
        
        rc = inst.SetAdcCaptureEnable(on);
        assert(rc == 0);
        
        rc = inst.ReadAdcCaptureStatus();
        
        rc = inst.ReadMultipleAdcFrames(chanIndex, frameNumber, numberOfFrames, netArray);
        assert(rc == 0);
        
        samples = uint16(netArray);
        
        plot(samples);
        axis([0 readLen 0 2^12]);
        yline(2047);
        
        pause(0.1);
        
        rc = inst.SetAdcCaptureEnable(off);
        assert(rc == 0);
        
    end   
    
    % Free the memory space that was allocated for ADC capture
    rc = inst.FreeAdcReservedSpace();
    assert(rc == 0);
    
    fprintf('ADC mem cleared\n');

% ---------------------------------------------------------------------
% End of the example
% ---------------------------------------------------------------------
    
%     Disable output
    res = inst.SendScpi(':OUTP OFF');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd of Example - %s\n', netStrToStr(res.RespStr));
    close all % Close all figures
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
%% Signal Processing Functions

%Create a sine wave
function sine = sine(cycles, phase, segLen, bits)
   
  verticalScale = ((2^bits))-1; % 2584 16 bits , 9082 8 bits
  time = -(segLen-1)/2:(segLen-1)/2;
  omega = 2 * pi() * cycles;
  rawSignal = sin(omega*time/segLen); 
  %rawSine = amp* cos(omega*time/segLen); 
 
  sine = ampScale(bits, rawSignal);
  
  %plot(sine);
  
end

function sine = addSinePulse(segment, starttime, dt, pulseDuration, freq, phase, bits)            
    rawSignal = segment;
    npoints = length(segment);
    starttimeIdx = roundCheck(starttime, dt)+1
    durationPts = roundCheck(pulseDuration, dt)
    
    if starttimeIdx > 0 && starttimeIdx + durationPts <= npoints
        rawSignal(starttimeIdx:starttimeIdx+durationPts) = sin(2*pi*freq*dt*(0:durationPts) + 2*pi*freq*starttime + phase);
        sine = ampScale(bits, rawSignal);
    else
        error('Pulse out of bounds');
    end
    figure(2)
    plot(sine)
end

function dacWav = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits)            

    t = 0:1/sampleRateDAC:rampTime;
    dacWave = chirp(t,fStart,rampTime,fStop);
    seglenTrunk = (floor(length(dacWave)/ 64))*64;
    dacWave = dacWave(1:seglenTrunk);
    dacWav = ampScale(bits, dacWave);

end

% Scale to FSD
function dacSignal = ampScale(bits, rawSignal)
 
  maxSig = max(rawSignal);
  verticalScale = ((2^bits)/2)-1;

  vertScaled = (rawSignal / maxSig) * verticalScale;
  dacSignal = uint8(vertScaled + verticalScale);
  %plot(dacSignal);

%   if bits > 8
%       dacSignal16 = [];
%       sigLen = length(dacSignal);
%       k=1;
%       for j = 1:2:sigLen*2;
%         dacSignal16(j) = bitand(dacSignal(k), 255);
%         dacSignal16(j+1) = bitshift(dacSignal(k),-8);
%         k = k + 1;
%       end
%       dacSignal = dacSignal16;
%   end
end

function rNb = roundCheck(nb, units)
    global debug
    rNb = round(nb/units);
    if ~isempty(debug) && debug > 1 && abs(rNb*units - nb) > units/100
       disp([inputname(1) ' = ' num2str(nb/units) ' rounded to ' num2str(rNb)]);
    end
end
