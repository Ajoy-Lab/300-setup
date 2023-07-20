%% Proteus P2584M Example
% Connect to P2584M instrument, and download waveform.

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

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
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    res = inst.SendScpi('*CLS');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':FREQ:RAST 2.5e9');
    assert(res.ErrCode == 0);            
    
    % ---------------------------------------------------------------------
    % Download sine waveform of 1024 points to segment 1 of channels 1 & 2
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);
    
    % Define segment 1 of 1024 points
    res = inst.SendScpi(':TRAC:DEF 1,1024');
    assert(res.ErrCode == 0);
    
    % build sinusoid wave of 1024 points
    low_dac_level = 0;
    high_dac_level = 65535;
    
    amp = (high_dac_level - low_dac_level) / 2.0;
    mid = (high_dac_level + low_dac_level) / 2.0;
    
    x = linspace(0, 2 * pi, 1024 + 1);
    wavdat = sin(x(1:1024)) * amp + mid;
    
    % convert to uint16
    wavdat = round(wavdat);
    wavdat = min(wavdat, high_dac_level);
    wavdat = max(wavdat, low_dac_level);    
    wavdat = uint16(wavdat);
    
    % get the bytes of the waveform data
    wavdat = typecast(wavdat, 'uint8');
    
    % select segmen 1 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 1');
    assert(res.ErrCode == 0); 
    
    % Download the binary data to segment 1
    res = inst.WriteBinaryData(':TRAC:DATA 0,#', wavdat);
    assert(res.ErrCode == 0);
    
    % ---------------------------------------------------------------------
    % Load waveform of 1024 points from file to segment 2 of channels 1 & 2
    % ---------------------------------------------------------------------
    
    % Define segment 1 of 1024 points
    res = inst.SendScpi(':TRAC:DEF 2,1024');
    assert(res.ErrCode == 0);
    
    % Get the path of the m file
    filepath = matlab.desktop.editor.getActiveFilename;
    % Extract the path of the directory
    [dirpath,fname,ext] = fileparts(filepath);
    
    % Build the path of the waveform file
    filepath = [ dirpath filesep 'triangle1024ptsDac16.wav' ];    
    % And convert it to byte array
    filepath = uint8(filepath);
    
    % Select segmen 2 as the the programmable segment
    res = inst.SendScpi(':TRAC:SEL 2');
    assert(res.ErrCode == 0); 
    
    % Load the binary data from the file to to segment 2
    res = inst.WriteBinaryData(':TRAC:FNAM 0,#', filepath);
    assert(res.ErrCode == 0);
    
    
    % ---------------------------------------------------------------------
    % Play segment 1 in channel 1
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);    
    
    res = inst.SendScpi(':SOUR:FUNC:SEG 1');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SOUR:VOLT 1.0');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);
    
    % ---------------------------------------------------------------------
    % Play segment 2 in channel 2
    % ---------------------------------------------------------------------
    
    res = inst.SendScpi('INST:CHAN 2');
    assert(res.ErrCode == 0);    
    
    res = inst.SendScpi(':SOUR:FUNC:SEG 2');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':SOUR:VOLT 1.0');
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi(':OUTP ON');
    assert(res.ErrCode == 0);
    
    % ---------------------------------------------------------------------
    % Create Task-Table
    % ---------------------------------------------------------------------
    
    taskTableLen = 2;
    taskTableRow = TaskInfo();
    rowBinarySize = taskTableRow.SerializedSize;
    tableBinDat = NET.createArray('System.Byte',taskTableLen * rowBinarySize);    
    
    % ----------
    % first row
    % ----------
    taskTableRow.SegNb = 1; 
    
    % TaskState is either
    %  TaskStateType.Single, 
    %  TaskStateType.StartSequence,
    %  TaskStateType.EndSequence or
    %  TaskStateType.InsideSequence
    taskTableRow.TaskState = TaskStateType.Single;
    
    taskTableRow.TaskLoopCount = 1;
    taskTableRow.SeqLoopCount = 1;
    
    % TaskIdleWaveform is either
    %  IdleWaveform.DC,
    %  IdleWaveform.FirstPoint or
    %  IdleWaveform.CurrentSeg
    taskTableRow.TaskIdleWaveform = IdleWaveform.DC;
    
    taskTableRow.TaskDcVal = 32768; % Mid-DAC (16 bits)
    
    % TaskEnableSignal is (currently) either
    %   EnableSignalType.None,
    %   EnableSignalType.ExternTrig1,
    %   EnableSignalType.ExternTrig2,
    %   EnableSignalType.InternTrig,
    %   EnableSignalType.CPU,
    %   EnableSignalType.FeedbackTrig or
    %   EnableSignalType.HwControl
    taskTableRow.TaskEnableSignal = EnableSignalType.None;
    
    % TaskAbortSignal is (currently) either
    %   AbortSignalType.None,
    %   AbortSignalType.ExternTrig1,
    %   AbortSignalType.ExternTrig2,
    %   AbortSignalType.InternTrig,
    %   AbortSignalType.CPU,
    %   AbortSignalType.FeedbackTrig or
    %   AbortSignalType.AnyExternTrig
    taskTableRow.TaskAbortSignal = AbortSignalType.None;
    
    % TaskAbortJump is either
    %   ReactMode.Eventually or
    %   ReactMode.Immediately
    taskTableRow.TaskAbortJump = ReactMode.Eventually;
    
    % TaskCondJumpDest is either
    %   TaskCondJump.Next1Task
    %   TaskCondJump.FeedbackTrigValue
    %   TaskCondJump.SwitchNext1Next2
    %   TaskCondJump.NextTaskSel
    %   TaskCondJump.NextScenario
    taskTableRow.TaskCondJumpDest = TaskCondJump.Next1Task;  
    
    taskTableRow.NextTask1 = 2;
    taskTableRow.NextTask2 = 0;
    
    taskTableRow.NextTaskDelay = 0;
    taskTableRow.TaskLoopTrigEnable = false;    
    
    % Serialize it at offset 0
    taskTableRow.Serialize(tableBinDat, 0);
    
    % ----------
    % second row
    % ----------
    
    taskTableRow.SegNb = 2;
    taskTableRow.TaskLoopCount = 4;
    taskTableRow.NextTask1 = 1;
    
    % Serialize it at offset rowBinarySize.
    % The offset of the n-th row is: (n-1)*rowBinarySize
    taskTableRow.Serialize(tableBinDat, int32(rowBinarySize));
    
    % ---------------------------------------------
    % Write the binary task-table data to channel 1
    % ---------------------------------------------
    
    res = inst.SendScpi('INST:CHAN 1');
    assert(res.ErrCode == 0);
    
    % optionaly wipe the task-table first
    res = inst.SendScpi('TASK:DEL:ALL');
    assert(res.ErrCode == 0);
    
    % Download binary-data
    res = inst.WriteBinaryData(':TASK:DATA 0,#', tableBinDat);
    assert(res.ErrCode == 0);
    
    % Set function mode = 'TASK'
    res = inst.SendScpi('FUNC:MODE TASK');
    assert(res.ErrCode == 0);
    
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
    catch ME
       rethrow(ME) 
    end
end


