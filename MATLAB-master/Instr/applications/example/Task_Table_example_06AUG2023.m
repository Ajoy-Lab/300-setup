% EXAMPLE FOR TASK TABLE CREATION AND DOWNLOAD
%=============================================
% This example calculates up to 4 different signals and download them into
% four consecutive segments in the target channel of the Proteus device to 
% be used in a sequence defined by the corresponding task table
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
% In this example, the Task Composer is used to create the table in the
% target Proteus. The Task Composer uses SCPI commands to define each field
% for each task table entry.

clc;

fprintf(1, 'INITIALIZING SETTINGS\n');

% Communication Parameters
remotePort = 2020;
localPort = 9090;

off = 0;
on = 1;
pfunc = ProteusFunctions;

dll_path = 'C:\\Windows\\System32\\TEPAdmin.dll';

cType = "DLL";  %"LAN" or "DLL"

paranoia_level = 2;

if cType == "LAN"
    try
        connStr = strcat('TCPIP::',connStr,'::5025::SOCKET');
        inst = TEProteusInst(connStr, paranoia_level);
        
        res = inst.Connect();
        assert (res == true);
    catch ME
        rethrow(ME)
    end   
else
    asm = NET.addAssembly(dll_path);

    import TaborElec.Proteus.CLI.*
    import TaborElec.Proteus.CLI.Admin.*
    import System.*
    
    admin = CProteusAdmin(@OnLoggerEvent);
    rc = admin.Open();
    assert(rc == 0);   
    
    try
        slotIds = admin.GetSlotIds();
        numSlots = length(size(slotIds));
        assert(numSlots > 0);
        
        % If there are multiple slots, let the user select one ..
        sId = slotIds(1);
        if numSlots > 1
            fprintf('\n%d slots were found\n', numSlots);
            for n = 1:numSlots
                sId = slotIds(n);
                slotInfo = admin.GetSlotInfo(sId);
                if ~slotInfo.IsSlotInUse
                    modelName = slotInfo.ModelName;
                    if slotInfo.IsDummySlot
                        fprintf(' * Slot Number:%d Model %s [Dummy Slot].\n', sId, modelName);
                    else
                        fprintf(' * Slot Number:%d Model %s.\n', sId, modelName);
                    end
                end
            end
            pause(0.1);
            choice = 8%input('Enter SlotId ');
            fprintf('\n');
            sId = uint32(choice);
        end
        
        % Connect to the selected instrument ..
        should_reset = true;
        inst = admin.OpenInstrument(sId, should_reset);
        instId = inst.InstrId;
        
    catch ME
        admin.Close();
        rethrow(ME) 
    end    
end

%% Create Administrator
% fprintf('\n');
% 
% res = inst.Connect();
% assert (res == true);

% Identify instrument using the standard IEEE-488.2 Command
% idnstr = inst.identifyModel();
% fprintf('\nConnected to: %s\n', idnstr);

% Reset AWG
inst.SendScpi('*CLS');
inst.SendScpi('*RST');

% Get options using the standard IEEE-488.2 Command
% optstr = inst.getOptions();

% Get Number of Channels
numOfChannels = 4
% This example is written to handle a maximum of 4 channels
if numOfChannels > 4
    numOfChannels = 4;
end

% Get maximum sample rate for target instrument
samplingRate = 1.0e9;

fprintf(1, 'Calculating WAVEFORMS\n');

% select channel in the range 1..4
channel = 1;
channel = mod(channel, numOfChannels);
% First segment to be used and total number of segments
segment = 1;
numOfSegments = 4;
% Num of cycles and period in the basic square waveform
minCycles = 1;
period = 1.0E-6;

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Set sampling rate for AWG to maximum.
inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);
% Get granularity
granul = 32;
% Get the default DAC resolution for the current settings.
dacRes = 16;

% Calculate basic square wave
myWfm = getSquareWfm(   samplingRate,... 
                        minCycles,...
                        period,...
                        granul);
                    
%Select Channel
inst.SendScpi(sprintf(':INST:CHAN %d', channel));
% DAC Mode set to 'DIRECT" (Default)
inst.SendScpi(':MODE DIR');
% All segments deleted
inst.SendScpi(':TRAC:DEL:ALL');

% Waveform Downloading
% ********************

for segNumber = 0:(numOfSegments - 1)
    
    fprintf(1, 'DOWNLOADING WAVEFORM FOR CH%d\n', channel);
    res = SendWfmToProteus(inst, channel, segment + segNumber, myWfm, dacRes);
    fprintf(1, 'WAVEFORM DOWNLOADED!\n');
        
    % The new waveform is calculated for the next channel
    if channel < numOfChannels
        % Integration, DC removal, and Normalization to the -1.0/+1.0 range
        myWfm = cumsum(myWfm);       
        myWfm = myWfm - mean(myWfm); 
        myWfm = myWfm / max(abs(myWfm));
    end
end

fprintf(1, 'CREATING SEQUENCE\n');

% Simple Sequence with N tasks
totalTasks = 4;
% The Task Composer is configured to handle a certain number of task
% entries
inst.SendScpi(sprintf(':TASK:COMP:LENG %d', totalTasks)); 

% Then, each task is defined
for taskNumber = 1:totalTasks
    % Task to be defined is selected
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d', taskNumber));
    % The type of task is defined. SINGle is the default so sending this 
    % command is not mandatory
    inst.SendScpi(':TASK:COMP:TYPE SING');
    % The action to take after completing the task is defined. NEXT is the
    % default so sending this command is not mandatory
    inst.SendScpi(':TASK:COMP:DEST NEXT');
    % Assigns segment for task in the sequence 1..numOfSegments
    param = mod(taskNumber - 1, numOfSegments) + 1;
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', param));    
    % Next Task is the following task in the table except for the last
    % task, which points to task #1
    param = mod(taskNumber, totalTasks) + 1;
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d', param));
    % Num of Loops = task #    
    param = taskNumber;     
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', param));
end

% The task table created with the Composer is written to the actual task
% table of teh selected channel
inst.SendScpi(':TASK:COMP:WRITE');
fprintf(1, 'SEQUENCE CREATED!\n');

fprintf(1, 'SETTING AWG OUTPUT\n');
% Select Task Mode for generation
inst.SendScpi(':FUNC:MODE TASK');
% Start in task #1 (#1 is the default)
inst.SendScpi(':FUNC:MODE:TASK 1')

% Output volatge set to MAX
inst.SendScpi(':SOUR:VOLT MAX');   
% Activate outpurt and start generation
inst.SendScpi(':OUTP ON');

% It is recommended to disconnect from instrument at the end
admin.CloseInstrument(inst.InstrId);
admin.Close();  
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
    instHandle.SendScpi(sprintf(':INST:CHAN %d', channel));    
    instHandle.SendScpi(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)));        
    % select segmen as the the programmable segment
    instHandle.SendScpi(sprintf(':TRAC:SEL %d', segment));

    % format Wfm
    myWfm = myQuantization(myWfm, dacRes);
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';
    
    if dacRes == 16
        myWfm = typecast(myWfm, 'uint8');
        instHandle.WriteBinaryData(prefix, myWfm);
    else
        myWfm = uint8(myWfm)
        instHandle.WriteBinaryData(prefix, myWfm);
    end   
    
    result = length(myWfm);
end  