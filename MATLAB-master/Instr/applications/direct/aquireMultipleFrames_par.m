clear;

%% ADC Example
% Capture Data with ADC

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('c:\windows\system32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

%% Create Administrator
% Create instance of the CProteusAdmin class, and open it.
% Note that only a single CProteusAdmin instance can be open, 
% and it must be kept open till the end of the session.

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
    should_reset = false;
    inst = admin.OpenInstrument(sId);
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
    
    %
    fprintf('Reset complete\n');
    
    % ---------------------------------------------------------------------
    % Operate the ADC with direct functions (tentative)
    % ---------------------------------------------------------------------
    
    
    sampleRate = 1000e6;
    %memoryAlloc = 6400000000;

    % ---------------------------------------------------------------------
    % ADC Config
    % ---------------------------------------------------------------------
    
    rc = inst.SetAdcDualChanMode(1);
    assert(rc == 0);
    
    rc = inst.SetAdcSamplingRate(sampleRate);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(0,800);
    assert(rc == 0);

    rc = inst.SetAdcCaptureTrigSource(2);
    assert(rc == 0);
    
    rc = inst.SetAdcSelfTrigThresh(0, 0.055);
    assert(rc == 0);
    

    numberOfPulses = 2717*4; % 1 second of pulses
    seconds = 2; % Total capture time
    readLen = 1056*34;
    offLen = 0;
    
    rc = inst.SetAdcAcquisitionEn(1,0);
    assert(rc == 0);
    
    rc = inst.SetAdcFramesLayout(numberOfPulses*seconds, readLen);
    assert(rc == 0);
    
    % ---------------------------------------------------------------------
    % Capture Parmeters
    % ---------------------------------------------------------------------
       
    rc = inst.SetAdcMultiFrameModeEn(1);
    assert(rc == 0);
     
    % ---------------------------------------------------------------------
    % ADC Capture
    % ---------------------------------------------------------------------
              
    rc = inst.SetAdcCaptureEnable(1);
    assert(rc == 0);
    choice = input('Wait for shuttle to complete then press [Enter] to continue ');
    fprintf(1, '\n');
    % Wait till the capture completes
    
    status = inst.ReadAdcCaptureStatus();
    for i = 1 : 250
        if status ~= 0
            break;
        end
        pause(0.01);
        status = inst.ReadAdcCaptureStatus();
    end
    
    %readSize = uint64(100000);
    %readOffset = uint64(100000);
    
    readSize = uint64(readLen);
    readOffset = uint64(offLen);
    
    
    chanIndex = 0;
    %netArray = NET.createArray('System.UInt16', 100000);
    
    netArray = NET.createArray('System.UInt16', readLen*numberOfPulses);
    
    %pulseAmp = [numberOfPulses*seconds];
    pulseAmp =[];
    pulseAmp_temp=zeros(1,numberOfPulses);
    padded_len=2^18;
    dF = sampleRate/padded_len;
    f = -sampleRate/2:dF:sampleRate/2-dF;
    
    [v,b1]=min(abs(f-20e6));
    [v,b2]=min(abs(f+20e6));
    
    eff_read=100:readLen-100;
    cyclesPoints = 50;
    
    for n = 1:seconds
        n
        fprintf('Start Read\n');
        rc = inst.ReadMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, netArray);
        assert(rc == 0);
        fprintf('Read Complete\n');

        fprintf('Start Processing\n');
        samples = double(netArray);
        tic
        
        for i = 1:numberOfPulses
            tic
%           clear pulse pulseMean rc X linFFT
            pulse = samples((readLen*(i-1))+1:readLen*(i));
            pulseMean = mean(pulse);
            pulseDC = pulse - pulseMean; % remove DC
            toc
            
            tic
            X = fftshift(fft(pulseDC));
            linFFT = (abs(X)/readLen);
            [amp, loc] = max(linFFT);
             pulseAmp_temp(i) = amp;
             toc
        end
        pulseAmp=[pulseAmp pulseAmp_temp];

        toc
        fprintf('End\n');
    end

    ivec=1:numberOfPulses*seconds;
    time_axis=92e-6.*ivec;
    figure(1);clf;
    plot(time_axis,pulseAmp);;
    set(gca,'ylim',[0 max(pulseAmp)*1.1]);    
    
    rc = inst.SetAdcCaptureEnable(0);
    assert(rc == 0);
    
    
    % Free the memory space that was allocated for ADC capture
    rc = inst.FreeAdcReservedSpace();
    assert(rc == 0);
    
    % ---------------------------------------------------------------------
    % End of the example
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd of Example - %s\n', netStrToStr(res.RespStr));
    
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
