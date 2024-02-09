% EXAMPLE FOR BINARY TASK TABLE CREATION AND DOWNLOAD
%====================================================
% This example calculates up to 4 different signals and download them into
% four consecutive segments in the target channel of the Proteus device to 
% be used in a sequence defined by the corresponding task table.
% 
% A task table is defined so all the four waveforms are sequenced in
% a continuous loop. In order to do so, each task must jump to the next
% excepto for the last one that must jump to task #1. Each task is looped
% (task #) number of times. This is the loop sequence:
%      
%  Task #1, Segm #1, 1 loop
%  Task #2, Segm #2, 2 loops
%  Task #3, Segm #3, 3 loops
%  Task #4, Segm #4, 4 loops
%
%  -->Task#1(1 loop)->Task#2(2 loop)->Task#3(3 loop)->Task#4(4 loop)->
% |                                                                   |
%  <------------------------------------------------------------------
%
% In this example, a task table entry object is defined, then the complete
% table is defined as a array of task table entry objects. For dwonload,
% the array of task table entries is converted to a pure sequence of bytes
% and downloaded as a binary block to the target Proteus device.

clc;

fprintf(1, 'INITIALIZING SETTINGS\n');

% Communication Parameters
connStr = '192.168.1.48'; % your IP here
paranoia_level = 2; % 0, 1 or 2

%% Create Administrator
inst = TEProteusInst(connStr, paranoia_level);
fprintf('\n');

res = inst.Connect();
assert (res == true);

% Identify instrument using the standard IEEE-488.2 Command
idnstr = inst.identifyModel();
fprintf('\nConnected to: %s\n', idnstr);

% Reset AWG
inst.SendCmd('*CLS');
inst.SendCmd('*RST');

% Get options using the standard IEEE-488.2 Command
optstr = inst.getOptions();

% Get Number of Channels
numOfChannels = inst.getNumOfChannels(idnstr);
% This example is written to handle a maximum of 4 channels
if numOfChannels > 4
    numOfChannels = 4;
end

% Get maximum sample rate for target instrument
samplingRate = inst.getMaxSamplingRate();

fprintf(1, 'Calculating WAVEFORMS\n');

channel = 1;
segment = 1;
numOfSegments = 4;

minCycles = 1;
period = 1.0E-6;

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Set sampling rate for AWG to maximum.
inst.SendCmd([':FREQ:RAST ' num2str(samplingRate)]);

% Get granularity
granul = inst.getGranularity(idnstr, optstr);

% Get the default DAC resolution for the current settings.
dacRes = inst.getDacResolution();

% Calculate basic square wave
myWfm = getSquareWfm(   samplingRate,... 
                        minCycles,...
                        period,...
                        granul);
                    
%Select Channel
inst.SendCmd(sprintf(':INST:CHAN %d', channel));
% DAC Mode set to 'DIRECT' (Default)
inst.SendCmd(':MODE DIR');
% All segments deleted
inst.SendCmd(':TRAC:DEL:ALL');

% Waveform Downloading
% *******************

for segNumber = 0:(numOfSegments - 1)  
    
    fprintf(1, 'DOWNLOADING WAVEFORM FOR SEGMENT #%d\n', segNumber);
    res = SendWfmToProteus(inst, channel, segment + segNumber, myWfm, dacRes);        
    % The new waveform is calculated for the next channel
    if channel < numOfChannels
        % Integration
        myWfm = cumsum(myWfm);
        % DC removal
        myWfm = myWfm - mean(myWfm);
        % Normalization to the -1.0/+1.0 range
        myWfm = myWfm / max(abs(myWfm));
    end
end

fprintf(1, 'WAVEFORMS DOWNLOADED!\n');

fprintf(1, 'CREATING SEQUENCE\n');

% Simple Sequence with N tasks
totalTasks = 4;
mySequence = CreateTaskTable(totalTasks);  

for taskNumber = 1:totalTasks
    % Assigns segment for task in the sequence 1..numOfSegments
    param = mod(taskNumber - 1, numOfSegments) + 1;
    mySequence(taskNumber) = SetValueInTask(mySequence(taskNumber),...
                        'segNb', param);  
    % Next Task is the next task entry except for the last one that points
    % to the first task.
    param = mod(taskNumber, totalTasks) + 1;
    mySequence(taskNumber) = SetValueInTask(mySequence(taskNumber),...
                        'nextTask1', param);    
    % Number of loops = task #
    param = taskNumber;
    mySequence(taskNumber) = SetValueInTask(mySequence(taskNumber),...
                            'taskLoopCount', param);
   
end

% Convert task table to binary format for download
binSequence = TaskTableToBin(mySequence);

% Select task mode
inst.SendCmd(':FUNC:MODE TASK');

fprintf(1, 'SEQUENCE CREATED!\n');
fprintf(1, 'DOWNLOADING SEQUENCE!\n');

% Task binary data is downloaded
prefix = ':TASK:DATA';
inst.SendBinaryData(prefix, binSequence, 'uint8');

fprintf(1, 'SEQUENCE CREATED!\n');


fprintf(1, 'SETTING AWG OUTPUT\n');
% Select Task Mode for generation
inst.SendCmd(':FUNC:MODE TASK');
% Start in task #1 (#1 is the default)
inst.SendCmd(':FUNC:MODE:TASK 1')

% Output volatge set to MAX
inst.SendCmd(':SOUR:VOLT MAX');   
% Activate output and start generation
inst.SendCmd(':OUTP ON');

% It is recommended to disconnect from instrument at the end
inst.Disconnect();    
clear inst;
clear;
fprintf(1, 'END\n');

function sqrWfm = getSquareWfm( samplingRate,... 
                                numCycles,...
                                period,...
                                granularity)
                            
    wfmLength = round(numCycles * period *samplingRate);
    wfmLength = round(wfmLength / granularity) * granularity;
    
    period = wfmLength / numCycles;    
    sqrWfm = 0:(wfmLength - 1);    
    sqrWfm = square(sqrWfm * 2 * pi / period);    
                            
end

function retval = myQuantization (myArray, dacRes)
  
  minLevel = 0;
  maxLevel = 2 ^ dacRes - 1;  
  numOfLevels = maxLevel - minLevel + 1;
  
  retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
  retval = retval + minLevel;
  
  retval(retval > maxLevel) = maxLevel;
  retval(retval < minLevel) = minLevel;

end

function result = SendWfmToProteus( instHandle,...
                                    channel,...
                                    segment,...
                                    myWfm,...
                                    dacRes)

    %Select Channel
    instHandle.SendCmd(sprintf(':INST:CHAN %d', channel));    
    instHandle.SendCmd(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)));        
    % select segmen as the the programmable segment
    instHandle.SendCmd(sprintf(':TRAC:SEL %d', segment));

    % format Wfm
    myWfm = myQuantization(myWfm, dacRes);
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';
    
    if dacRes == 16
        instHandle.SendBinaryData(prefix, myWfm, 'uint16');
    else
        instHandle.SendBinaryData(prefix, myWfm, 'uint8');
    end   
    
    result = length(myWfm);
end  

%********************************************************
%           FUNCTIONS TO HANDLE BINARY TASK TABLES
%********************************************************

function outTaskTable = CreateTaskTable(n)   
    % outTaskTable - Create a task table with a number of entries
    %
    % Synopsis
    %   outTaskTable = CreateTaskTable(n)
    %
    % Description
    %   It creates a table of task table entry structs. The number of 
    %   entries is supplied as a paramter and the returned entries are
    %   initialized to the default values.
    %
    % Inputs ([]s are optional)
    %   (scalar) n   number of entries
    %
    % Outputs ([]s are optional)
    %   (task object) outTaskTablet   Array of Task Table structs
        
    outTaskTable = GetDefaultTask();
    outTaskTable = repmat(outTaskTable, 1, n);
end

function taskEntry = GetTaskEntry(  segN, ...
                                    nxTask1, ... 
                                    nxTask2, ...
                                    taskLoop, ...
                                    seqLoop, ...
                                    nextDelay, ...
                                    dcVal, ...
                                    idlWvf, ...
                                    enabSig, ...
                                    abortSig, ...
                                    condJump, ...
                                    jumpType,...
                                    state, ...
                                    loopTrig, ...
                                    adcTrig) 
                                
    % GetTaskEntry - Set all the values for a task struct
    %
    % Synopsis
    %   taskEntry = GetTaskEntry(....
    %
    % Description
    %   It sets all the fields of a task struct with a single function call
    %
    % Inputs ([]s are optional)
    %   See description in the function body
    %
    % Outputs ([]s are optional)
    %   (task object) taskEntry   Task struct with the values as specified
                                
    % The segment number
    taskEntry.segNb =               uint32(segN);
    % The Next-Task for Trigger 1 (zero for end)
    taskEntry.nextTask1 =           uint32(nxTask1);
    % The Next-Task for Trigger 2 (zero for end)
    taskEntry.nextTask2 =           uint32(nxTask2);
    % The task loop count (0:2^20-1)
    taskEntry.taskLoopCount =       uint32(taskLoop);
    % The sequence loop count (0:2^20-1).
    taskEntry.seqLoopCount =        uint32(seqLoop);
    % The delay in clocks before executing the next task.
    taskEntry.nextTaskDelay =       uint16(nextDelay);
    % The DAC value of the idle task DC waveform.  
    taskEntry.taskDcVal =           uint16(dcVal);
    % The behavior during idle-time
    % (0: DC, 1: First Point, 2: Current Segment(?))
    taskEntry.taskIdlWvf =          uint8(idlWvf);
    % The enabling signal
    % (0:None, 1:ExternTrig1, 2:ExternTrig2, 3:InternTrig, 
    %  4:CPU, 5:FeedbackTrig, 6:HW-Ctrl(?))
    taskEntry.taskEnableSig =       uint8(enabSig);
    % The aborting signal
    % (0:None, 1:ExternTrig1, 2:ExternTrig2, 3:InternTrig, 
    %  4:CPU, 5:FeedbackTrig, 6:Any)
    taskEntry.taskAbortSig =        uint8(abortSig);
    % How to decide where to jump
    % 0: Next1, 1: By FBTrig-Value, 2: ExtTrig[1/2]->Next[1/2],
    % 3: NextTaskSel(?), 4: Next Scenario)
    taskEntry.taskCondJumpSel =     uint8(condJump);
    % Task abort jump type
    % (0:Eventually, 1:Immediate)
    taskEntry.taskAbortJumpType =   uint8(jumpType);
    % The task state
    % (0:Single,1:First of sequence, 2:Last of sequence, 3:Inside Sequence)
    taskEntry.taskState =           uint8(state);
    % Enable (1) or disable (0) waiting for trigger on looping.
    taskEntry.taskLoopTrigEn =      uint8(loopTrig);
    % If it's non-zero, gen ADC trigger at the end of the current task.
    taskEntry.genAdcTrigger =       uint8(adcTrig);
end

function taskEntry = GetDefaultTask() 

    taskEntry = GetTaskEntry(   1, ... %Segment Number                               
                                1, ... %Next Task 1
                                1, ... %Next Task 2
                                1, ... %Task Loop
                                1, ... %Segment Loop
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0, ...
                                0);
end

function taskEntry = SetValueInTask(inpuTask, fieldName, value) 
                                
    taskEntry = inpuTask;
    fieldName = upper(fieldName);
    % This function allows for filed updating without taking care of
    % numeric types
    if strcmp(fieldName, upper('segNb'))
        taskEntry.segNb =               uint32(value);
    elseif strcmp(fieldName, upper('nextTask1'))
        taskEntry.nextTask1 =           uint32(value);
    elseif strcmp(fieldName, upper('nextTask2'))
        taskEntry.nextTask2 =           uint32(value);
    elseif strcmp(fieldName, upper('taskLoopCount'))
        taskEntry.taskLoopCount =       uint32(value);
    elseif strcmp(fieldName, upper('seqLoopCount'))
        taskEntry.seqLoopCount =        uint32(value);
    elseif strcmp(fieldName, upper('nextTaskDelay'))
        taskEntry.nextTaskDelay =       uint16(value);
    elseif strcmp(fieldName, upper('taskDcVal'))
        taskEntry.taskDcVal =           uint16(value);
    elseif strcmp(fieldName, upper('taskIdlWvf'))
        taskEntry.taskIdlWvf =          uint8(value);
    elseif strcmp(fieldName, upper('taskEnableSig'))
        taskEntry.taskEnableSig =       uint8(value);
    elseif strcmp(fieldName, upper('taskAbortSig'))
        taskEntry.taskAbortSig =        uint8(value);
    elseif strcmp(fieldName, upper('taskCondJumpSel'))
        taskEntry.taskCondJumpSel =     uint8(value);
    elseif strcmp(fieldName, upper('taskAbortJumpType'))
        taskEntry.taskAbortJumpType =   uint8(value);
    elseif strcmp(fieldName, upper('taskState'))
        taskEntry.taskState =           uint8(value);
    elseif strcmp(fieldName, upper('taskLoopTrigEn'))
        taskEntry.taskLoopTrigEn =      uint8(value);
    elseif strcmp(fieldName, upper('genAdcTrigger'))
        taskEntry.genAdcTrigger =       uint8(value);
    end
end

function taskTableBin = TaskTableToBin(taskEntry) 
                                
    taskTableBin = [];
    %Binary data is the concatenation of all the fields in all the entries
    %in the table formatted as bytes or 'uint8'
    for i = 1:length(taskEntry)
        val = taskEntry(i).segNb;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).nextTask1;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).nextTask2;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).taskLoopCount;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).seqLoopCount;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).nextTaskDelay;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).taskDcVal;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).taskIdlWvf;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).taskEnableSig;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];
        val = taskEntry(i).taskAbortSig;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];  
        val = taskEntry(i).taskCondJumpSel;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).taskAbortJumpType;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).taskState;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];    
        val = taskEntry(i).taskLoopTrigEn;
        taskTableBin = [taskTableBin typecast(val, 'uint8')]; 
        val = taskEntry(i).genAdcTrigger;
        taskTableBin = [taskTableBin typecast(val, 'uint8')];         
    end
end