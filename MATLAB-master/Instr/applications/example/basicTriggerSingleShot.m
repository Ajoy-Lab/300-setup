
%% ADC Example
% Capture Data with ADC

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
    
    fprintf('Reset complete\n');
    
    % ---------------------------------------------------------------------
    % Download sine waveform of 1024 points to segment 1 of channel 1
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);
    
    % Define segment 1 of 1024 points
    res = inst.SendScpi(':TRAC:DEF 1,1024');
    assert(res.ErrCode == 0);
    
    % build sinusoid wave of 1024 points
    low_dac_level = 0;
    high_dac_level = 255;
    
    amp = (high_dac_level - low_dac_level) / 2.0;
    mid = (high_dac_level + low_dac_level) / 2.0;
    
    x = linspace(0, 2 * pi, 1024 + 1);
    wavdat = 0.5*(sin(x(1:1024)) * amp + mid);
    
    % convert to uint8
    wavdat = round(wavdat);
    wavdat = min(wavdat, high_dac_level);
    wavdat = max(wavdat, low_dac_level);    
    wavdat = uint8(wavdat);
    
    % get the bytes of the waveform data
    wavdat = typecast(wavdat, 'uint8');
    
    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 
    
    % Download the binary data to segment 1
    res = inst.WriteBinaryData(':TRAC:DATA 0,#', wavdat);
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
    
    
    sampleRate = 800e6;
    %memoryAlloc = 6400000000
    memoryAlloc = 640000
    % Turn on ADC dual-channels mode (state = 1)
    rc = inst.SetAdcDualChanMode(1);
    assert(rc == 0);
    
    % Set ADC sampling  
    % which is he minimum in case of dual channels mode
    rc = inst.SetAdcSamplingRate(sampleRate);
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
    
    externalTrigger = false;
    
    if externalTrigger == true
        rc = inst.SetAdcExternTrigMode(1);
        assert(rc == 0);
        
        rc = inst.SetAdcExternTrigEn(1);
        assert(rc == 0);
        
        rc = inst.SetAdcExternTrigPulsDetWidth(0);
        assert(rc == 0);
        
        rc = inst.SetAdcExternTrigNegativePolarity(0);
        assert(rc == 0);
        
        choice = input('Generate external-trigger and press [Enter] to continue ');
        fprintf(1, '\n');
    else
        rc = inst.SetAdcExternTrigMode(0);
        assert(rc == 0);
        % Generate software-trigger for ADC capturing
        rc = inst.GenerateAdcSoftTrig();
        assert(rc == 0);        
    end
        
    
    % Wait till the capture completes
    status = inst.ReadAdcCaptureDoneStatus();
    for i = 1 : 250
        if status ~= 0
            break;
        end
        pause(0.01);
        status = inst.ReadAdcCaptureDoneStatus();
    end
    
    
    fprintf('Start Read\n');
    % Read x samples from offset y at the captured data of 
    % channel 1 to uint16 array and plot it
    
    readLen = 2000;
    offLen = 0 
    
    %readSize = uint64(100000);
    %readOffset = uint64(100000);
    
    readSize = uint64(readLen);
    readOffset = uint64(offLen);
    
    
    chanIndex = 0;
    %netArray = NET.createArray('System.UInt16', 100000);
    netArray = NET.createArray('System.UInt16', readLen);
    fprintf('End Read\n');
    
    rc = inst.ReadAdcChanData(chanIndex, readSize, readOffset,netArray);
    assert(rc == 0);
    
    samples = uint16(netArray);
    plot(samples);
    
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
