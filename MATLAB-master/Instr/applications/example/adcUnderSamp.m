%% ADC Example
% Capture Data with ADC

addpath '..\..\..\Utils';
addpath '..\..\..\Functions';

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Program Files\Tabor Electronics\ProteusGui\TEPAdmin.dll');

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
    % Download sine waveform of 1024 points to segment 1 of channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);
    
    sampleRateDAC = 9E9;
    segLen = 4096;
    bits = 8;
    amplitude = 1;
    cycles = 1092;

    % Define segment 1 
    res = inst.SendScpi(':TRAC:DEF 1,4096');
    assert(res.ErrCode == 0);
    
    dacSignal = amplitude * sine(sampleRateDAC, cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)

    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 

    % Download the binary data to segment 1
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
    % Operate the ADC with direct functions (tentative)
    % ---------------------------------------------------------------------


    sampleRateADC = 2e9;
    %memoryAlloc = 6400000000
    memoryAlloc = 6400

    readLen = 2000;
    offLen = 0 

    %readSize = uint64(100000);
    %readOffset = uint64(100000);

    readSize = uint64(readLen);
    readOffset = uint64(offLen);

    chanIndex = 0;
    %netArray = NET.createArray('System.UInt16', 100000);
    netArray = NET.createArray('System.UInt16', readLen);


    % Turn on ADC dual-channels mode (state = 1)
    rc = inst.SetAdcDualChanMode(1);
    assert(rc == 0);

    % Set ADC sampling  
    % which is he minimum in case of dual channels mode
    rc = inst.SetAdcSamplingRate(sampleRateADC);
    assert(rc == 0);

    % Get the granularity of the capture-size and capture-offset in samples
    % per channel (both depend on the state of the dual-channels mode)    
    sizeGranularity = inst.GetAdcCaptureSizeGranularity();
    offsetGranularity = inst.GetAdcCaptureOffsetGranularity();

    % Allocate memory-space for ADC capture.
    % This space is taken from the space of the waveform-segments
    % therefore, one should free the ADC space when it is not needed.
    %
    % Note that the size of the allocated-space (in samples per channel) 
    % should bean integral multiply of the capture-granularity, and if it
    % is not then it will be automatically increased to the nearest
    % multiply of the capture-granularity
    rc = inst.AllocAdcReservedSpace(memoryAlloc);
    assert(rc == 0);

    % Optionally query the size of the space that was allocated
    allocSpace = inst.GetAdcReservedSpaceSize();

    % Set the capture size and offset inside the reserved space.
    % Both expressed in samples per channel, and should be multiplies
    % of the corresponding granularities.
    %
    % Note that if the capture-size is not a multiply of the
    % size-granularity, then it will be increased to the nearest multiply,
    % but if the offset is not a multiply of the offset granularity, then 
    % the function fails and returns error.    
    rc = inst.SetAdcCaptureOffset(0);
    assert(rc == 0);

    rc = inst.SetAdcCaptureSize(memoryAlloc);
    assert(rc == 0);

    % Optionally query the capture size
    captureSize = inst.GetAdcCaptureSize();

    fprintf('Capture size set\n');
    h = figure(1);
    forever = 1; 
    while forever 
        drawnow
        isKeyPressed = ~isempty(get(h,'CurrentCharacter'));
        if isKeyPressed
            break
        end
        rc = inst.SetAdcExternTrigMode(0);
        assert(rc == 0);
        % Generate software-trigger for ADC capturing
        rc = inst.GenerateAdcSoftTrig();
        assert(rc == 0);        

        % Wait till the capture completes
        status = inst.ReadAdcCaptureDoneStatus();
        for i = 1 : 2500
            if status ~= 0
                break;
            end
            pause(0.01);
            status = inst.ReadAdcCaptureDoneStatus();
        end

        rc = inst.ReadAdcChanData(chanIndex, readSize, readOffset,netArray);
        assert(rc == 0);

        samples = int16(netArray);
        
       
        maxSig = max(samples);
        minSig = min(samples);
        meanSig = mean(samples);
        signalRangePkPk = (minSig*-1)+maxSig;

        dataReadTimeDC=samples-meanSig;

        X = fftshift(fft(dataReadTimeDC));
        N=readLen;
        dF = sampleRateADC/N;                      % hertz
        f = -sampleRateADC/2:dF:sampleRateADC/2-dF;

        tiledlayout(1,3);
        nexttile;
        plot(samples);
        title('Waveform');
        xlabel('Points')
        nexttile;
        plot(dataReadTimeDC);
        title('Waveform DC Removed');
        xlabel('Points');
        nexttile;
        plot(f,abs(X)/N);
        maxpeak = max(abs(X)/N);
        title('Spectrum');
        xlabel('Hz');

        pause(0.05);
    end    

    % Read The 1G samples to binary file
    %chanIndex = 0;
    %readSize = uint64(1000000000);
    %readOffset = uint64(0);
    %fileWrOffset = uint64(0);
    %filePath = 'C:\Tabor\Customers\Berkeley\Customeradc_capture.wav';

    %rc = inst.ReadAdcChanData(chanIndex, readSize, readOffset, filePath, fileWrOffset);
    %assert(rc == 0);


    % Free the memory space that was allocated for ADC capture
    %rc = inst.FreeAdcReservedSpace();
    %assert(rc == 0);

% ---------------------------------------------------------------------
% End of the example
% ---------------------------------------------------------------------
    
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

