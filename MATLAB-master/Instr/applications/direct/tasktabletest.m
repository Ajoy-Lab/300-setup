%% Initializing AWG Receiving & Digitizing Functionality
clear;clc;

pfunc = ProteusFunctions;

fprintf(1, 'INITIALIZING SETTINGS\n');

dll_path = 'C:\\Windows\\System32\\TEPAdmin.dll';
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
    sId = 8;
    % Connect to the selected instrument ..
    should_reset = true;
    inst = admin.OpenInstrument(sId, should_reset);
    instId = inst.InstrId;

catch ME
    admin.Close();
    rethrow(ME) 
end

% ---------------------------------------------------------------------
% Setup instrument
% ---------------------------------------------------------------------

res = inst.SendScpi('*IDN?'); %asking for instrument IP address, model/serial #
assert(res.ErrCode == 0);
fprintf(1, '\nConnected to ''%s''\n', pfunc.netStrToStr(res.RespStr));


% SET UP AWG
fprintf(1, 'SETTING AWG\n');

awgSRate = 1E9
awgChan = 2

%Select Channel

inst.SendScpi("*CLS; *RST")
inst.SendScpi(sprintf(':INST:CHAN %d', awgChan));
inst.SendScpi([':FREQ:RAST ' num2str(awgSRate)]);
inst.SendScpi(":INIT:CONT OFF");

% sequence 1 creation
chirps{1}.cycles = 16;
chirps{1}.phase = 0;
chirps{1}.segLen = 20480;
chirps{1}.bits = 8;
chirps{1}.dacSignal = sine(chirps{1}.cycles, chirps{1}.phase, chirps{1}.segLen, chirps{1}.bits);
chirps{1}.segm = 1;

% sequence 2 creation
chirps{2}.cycles = 2;
chirps{2}.phase = 0;
chirps{2}.segLen = 20480;
chirps{2}.bits = 8;
chirps{2}.dacSignal = sine(chirps{2}.cycles, chirps{2}.phase, chirps{2}.segLen, chirps{2}.bits);
chirps{2}.segm = 2;

pol_times = [8e-3,0.4e-3, 0.6e-3, 1e-3];

task_list = build_tasktable(inst,pol_times,chirps,'first sign','-');

% % Play seg 1
res = inst.SendScpi(':SOUR:FUNC:MODE:SEGM 1');
assert(res.ErrCode == 0);

res = inst.SendScpi(':SOUR:VOLT MAX');
assert(res.ErrCode == 0);

res = inst.SendScpi(':OUTP ON');
assert(res.ErrCode == 0);
% 
% % trig settings for ARB mode
inst.SendScpi(':TRIG:ACTIVE:SEL TRG2'); 
inst.SendScpi(':TRIG:LEV 1.0');
inst.SendScpi(':TRIG:ACTIVE:STAT ON');

% The Task Composer is configured to handle a certain number of task entries
% inst.SendScpi(':TASK:COMP:LENG 2'); 

% create the task table
create_task_table(inst,task_list, awgSRate);

% write task table
inst.SendScpi(':TASK:COMP:WRITE');
fprintf(1, 'SEQUENCE CREATED!\n');

fprintf(1, 'SETTING AWG OUTPUT\n');

inst.SendScpi(':SOUR:FUNC:MODE TASK');

res = inst.SendScpi(':OUTP ON');
assert(res.ErrCode == 0);


% ---------------------------------------------------------------------
% ADC Config
% ---------------------------------------------------------------------

sampleRate = 1000e6;
adcDualChanMode = 1;
fullScaleMilliVolts =1000;
adcChanInd = 0; % ADC Channel 1
trigSource = 1; % 1 = external-trigger
numberOfPulses = 10;
capture_first = 1;
capture_count = numberOfPulses;
loops = 1;
tacq = 32; % acquisition time as in pulse sequence
readLen = round2((tacq+2)*1e-6/1e-9,96)-96; % must be divisible by 96,
%netArray = NET.createArray('System.UInt16', readLen*numberOfPulses); %total array -- all memory needed
%number of points in a given small chunk

% Setup the digitizer 
res = inst.SendScpi(':DIG:MODE DUAL');

% Enable capturing data from channel 1
res = inst.SendScpi(':DIG:CHAN:SEL 1');

cmd = [':DIG:FREQ ' num2str(sampleRate)];
res = inst.SendScpi(cmd);

cmd = ':DIG:FREQ?';
res = inst.SendScpi(cmd);
fprintf(pfunc.netStrToStr(res.RespStr));

% Enable acquisition in the digitizer's channels  
res = inst.SendScpi(':DIG:CHAN:STATE ENAB');

% Setup frames layout   
cmd = [':DIG:ACQuire:FRAM:DEF ' num2str(numberOfPulses) ',' num2str(readLen)];
res = inst.SendScpi(cmd);

% Select internal or external as start-capturing trigger:
%res = inst.SendScpi(':DIG:TRIG:SOURCE CH1');
res = inst.SendScpi(':DIG:TRIG:SOURCE TASK1');

res = inst.SendScpi(':DIG:CHAN:SEL 1');

% Select which frames are filled with captured data 
%(all frames in this example)
inst.SendScpi(':DIG:ACQ:FRAM:CAPT:ALL');

% Delete all wafm memory
inst.SendScpi(':DIG:ACQ:ZERO:ALL');

resp = inst.SendScpi(':DIG:DATA:FORM?');
resp = strtrim(pfunc.netStrToStr(resp.RespStr));

inst.SendScpi(':DIG:INIT OFF'); 
    
res = inst.SendScpi(':DIG:CHAN:SEL 1');

inst.SendScpi(':DIG:INIT ON');  

for n = 1:150
    resp = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
    resp = strtrim(pfunc.netStrToStr(resp.RespStr));
    items = split(resp, ',');
    items = str2double(items);
    if length(items) >= 3 && items(2) == 1
        break
    end
    if mod(n, 10) == 0                
        fprintf('%d. %s Time:\n', fix(n / 10), resp);
%                 toc                
    end
    pause(0.1);
end
res = inst.SendScpi(':DIG:INIT OFF');
pause(1);

% After this point aquisition is complete, next stage is to readout data

% Define what we want to read(frames data, frame-header, or both).
res = inst.SendScpi(':DIG:DATA:TYPE FRAM');

% Choose which frames to read (all in this example)
res = inst.SendScpi(':DIG:DATA:SEL ALL');

% Read binary block
resp = inst.SendScpi(':DIG:DATA:SIZE?');
resp = strtrim(pfunc.netStrToStr(resp.RespStr));
num_bytes = str2double(resp);
% 
% Read the data that was captured by channel 1:
inst.SendScpi(':DIG:CHAN:SEL 1');

% because read format is UINT16 we divide byte number by 2
wavlen = floor(num_bytes / 2);

% allocate NET array
netArray = NET.createArray('System.UInt16', wavlen);
% read the captured frame
%tic
res = inst.ReadMultipleAdcFrames(0, capture_first, numberOfPulses, netArray);
%toc
assert(res == 0);

% cast to matlab vector
samples = uint16(netArray);

% deallocate the NET array
delete(netArray);

plot(samples);

% res = inst.SendScpi(':DIG:ACQ:ZERO:ALL 0');


% It is recommended to disconnect from instrument at the end
rc = admin.CloseInstrument(instId);
assert(rc==0);
fprintf('Instrument closed\n');

% Close the administrator at the end ..
admin.Close();

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

function u = create_task_table(inst,task_list, awgSRate)
limit_cycles = 1.04e6;
task_list = arrange_tasktable(task_list,awgSRate,limit_cycles);
    for i = 1:length(task_list)
        len = task_list{i}.len;
        dacSignal = task_list{i}.dacSignal;
        t = task_list{i}.time;
        segm = task_list{i}.segm;
        segm_period = len/awgSRate;
        cycles = round(t/segm_period);
%         cycles = 1;
        real_time = cycles*segm_period;
        inst.SendScpi(sprintf(':TASK:COMP:SEL %d',i));
        inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',segm));
        inst.SendScpi(':TASK:COMP:TYPE SING');
        inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',cycles));
        if i==1 % define trigger
            inst.SendScpi(':TASK:COMP:ENAB TRG2');
        end
        if i==length(task_list)
        else
            inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',i+1));
        end
    end
    u = 0;
end

function new_task_list = arrange_tasktable(task_list,awgSRate,limit_cycles)
num_task = 1;
     for i = 1:length(task_list)
            len = task_list{i}.len;
            dacSignal = task_list{i}.dacSignal;
            t = task_list{i}.time;
            segm = task_list{i}.segm;
            segm_period = len/awgSRate;
            cycles = round(t/segm_period);
            real_time = cycles*segm_period;
            if cycles >= limit_cycles
                n = fix(cycles/limit_cycles);
                for j = 1:n
                    new_task_list{num_task}.len = len;
                    new_task_list{num_task}.dacSignal = dacSignal;
                    new_task_list{num_task}.segm = segm;
                    new_task_list{num_task}.cycles = limit_cycles;
                    new_task_list{num_task}.time = limit_cycles*segm_period;
                    num_task = num_task +1;
                end
            end
            rem_cycles = rem(cycles, limit_cycles);
            if rem_cycles ~= 0
                 new_task_list{num_task}.len = len;
                 new_task_list{num_task}.dacSignal = dacSignal;
                 new_task_list{num_task}.segm = segm;
                 new_task_list{num_task}.cycles = rem_cycles;
                 new_task_list{num_task}.time = rem_cycles*segm_period;
                 num_task = num_task +1;
            end
     end
end

function task_list = build_tasktable(inst,pol_times,chirps,varargin)
% function that create the task liste before creating the task table
% input: - pol_times : list that contains the different DNP time for each
% steps. 
%          - chirps: 2 dimension structure that contains the information
%          for each chirp, with positive and negative DNP respectively
%
% options : -  'first sign' could be '+' or '-' to stqrt with positive or
% negative DNP. positive DNP is the default value

    chirp_num = 1;
    if isempty(varargin) == 0
        for i = length(varargin)/2
            if strcmp(varargin(i),'first sign') == 1
                if strcmp(varargin(i+1),'+') ==1
                    chirp_num = 1;
                elseif strcmp(varargin(i+1),'-')
                    chirp_num = 2;
                end
            end
        end
    end
    
    for i = 1:2
             % Define segment i 
        res = inst.SendScpi([':TRAC:DEF ',num2str(chirps{i}.segm),',',num2str(length(chirps{i}.dacSignal))]); % define memory location 1 of length dacSignal
        assert(res.ErrCode == 0);

        % select segment as the the programmable segment
        res = inst.SendScpi(sprintf(':TRAC:SEL %d',chirps{i}.segm));
        assert(res.ErrCode == 0); 

        % Download the binary data to segment
        res = inst.WriteBinaryData(':TRAC:DATA 0,', chirps{i}.dacSignal);
        assert(res.ErrCode == 0);
       
    end
    
    for i = 1:length(pol_times)
        task_list{i}.len = chirps{chirp_num}.segLen;
        task_list{i}.segm = chirps{chirp_num}.segm;
        task_list{i}.dacSignal = chirps{chirp_num}.dacSignal;
        task_list{i}.time = pol_times(i);
        if chirp_num == 1 
            chirp_num = 2;
        else
            chirp_num = 1;
        end
    end
end