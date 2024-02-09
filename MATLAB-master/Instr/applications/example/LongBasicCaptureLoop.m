clear;

%% ADC Example
% Capture Data with ADC

addpath 'C:\matlab-codes\MATLAB-master\MATLAB-master\Utils';
addpath 'C:\matlab-codes\MATLAB-master\MATLAB-master\Functions';

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Users\Simon\Documents\Tabor Electronics\ProteusAwg-10204-Customer\TEPAdmin.dll');

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
    
    res = inst.SendScpi(':FREQ:RAST 1E9'); % set sample clock
    assert(res.ErrCode == 0);    
    
    fprintf('Reset complete\n');

    % ---------------------------------------------------------------------
    % Download sine waveform of 1024 points to segment 1 of channel 1
    % ---------------------------------------------------------------------

    res = inst.SendScpi('INST:CHAN 1'); % select channel 1
    assert(res.ErrCode == 0);
    
    segLen = 2048;
    bits = 8;
    amplitude = 1;
    %cycles = 1092;
    cycles = 40;

    % Define segment 1 
    res = inst.SendScpi(':TRAC:DEF 1,2048'); % define memory location 1 as 2048 in length
    assert(res.ErrCode == 0);
    
    dacSignal = amplitude * sine(cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)

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

    res = inst.SendScpi(':SOUR:VOLT 0.4');
    assert(res.ErrCode == 0);

    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);

    fprintf('Waveform generated and playing\n');

% ---------------------------------------------------------------------
    % Operate the ADC with direct functions (tentative)
    % ---------------------------------------------------------------------

    sampleRateADC = 1e9;
    %memoryAlloc = 6400000000

    readLen = 33888; % integer div of 96
    offLen = 0 
    chanIndex = 0; %CH 1 of Digitizer
    off = 0;
    on = 1;
    numberOfFrames = 249980/4;
    frameNumber = 1;
    netArray = NET.createArray('System.UInt16', readLen);

    % Turn on ADC dual-channels mode (state = 1)
    rc = inst.SetAdcDualChanMode(on);
    assert(rc == 0);

    rc = inst.SetAdcSamplingRate(sampleRateADC);
    assert(rc == 0);

%     rc = inst.SetAdcMultiFrameModeEn(on);
%     assert(rc == 0);
    
    rc = inst.SetAdcDualChanMode(on);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(chanIndex,1000); % mv - 500, 800, 1000
    assert(rc == 0);

    rc = inst.SetAdcFramesLayout(numberOfFrames, readLen);
    assert(rc == 0);
    fprintf('Capture size set\n');
 
%     rc = inst.SetAdcCaptureTrigSource(chanIndex, 2); % 1 means external trigger, 2 means self trigger
%     assert(rc == 0);
    
%     rc = inst.SetAdcSelfTrigThresh(chanIndex, 0.05);
%     assert(rc == 0);
    

    trigSource = 1; % External
%     trigSource = 2; % Internal
 
    
    rc = inst.SetAdcCaptureTrigSource(chanIndex, trigSource); % 1 means external trigger, 2 means self trigger
    assert(rc == 0);
   
    rc = inst.SetAdcAcquisitionEn(on,off);
    assert(rc == 0);
    
    rc = inst.SetAdcExternTrigParams(0,1,0,1,0,0);
    assert(rc == 0); 

    rc = inst.SetAdcExternTrigThresh(0,0.5);
    assert(rc==0);

%     rc = inst.SetAdcExternTrigThresh(1,2.0);
%     assert(rc==0);
    


    fprintf('Trigger set\n');
   
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
        rc = inst.ReadMultipleAdcFrames(chanIndex, frameNumber, numberOfFrames, netArray);
        assert(rc == 0);
        
        samples = uint16(netArray);
        
%         ivec=1:readLen;
%         time_axis=ivec./sampleRateADC;
        plot(samples);
        
        pause(0.1);
    end   
    
    rc = inst.SetAdcCaptureEnable(off);
    assert(rc == 0);

    fprintf('Aquisition Disabled\n');
    
    
    % Free the memory space that was allocated for ADC capture
    rc = inst.FreeAdcReservedSpace();
    assert(rc == 0);
    
    fprintf('ADC mem cleared\n');

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

