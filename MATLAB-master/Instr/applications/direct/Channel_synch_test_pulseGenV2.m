% This version runs a Pulsed Spin-locking sequence
%% Clear everything
clear;
close;

%% Set defaults Vars
savealldata=false;
savesinglechunk=false;
savemultwind=false;

sampleRate = 2700e6;
global sampleRateDAC
sampleRateDAC = 9e9;
global inst
global interp
global pulseDict
global blockDict
pulseDict = containers.Map;
blockDict = containers.Map;
interp = 4;
adcDualChanMode = 2;
% fullScaleMilliVolts =1000;
trigSource = 1; % 1 = external-trigger
dacChanInd = 3;
adcChanInd = 2;
measurementTimeSeconds = 7; %Integer
%delay = 0.0000178; % dead time
%delay = 0.0000348; % dead time
%delay=0.00000148;
%delay = 0.0000108; % dead time
%delay = 0.0000028; % dead time
delay = 0.0000038; % dead time
global bits
bits = 16;


% remoteAddr = '192.168.1.2'; % old computer
remoteAddr = '192.168.10.5'; % new computer
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
    
% ---------------------------------------------------------------------
% Setup instrument
% ---------------------------------------------------------------------
%% case 1: initialize Tabor 
res = inst.SendScpi('*IDN?');
assert(res.ErrCode == 0);
fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));

pause(0.01);

res = inst.SendScpi('*CLS'); % clear
assert(res.ErrCode == 0);

res = inst.SendScpi('*RST'); % reset
assert(res.ErrCode == 0);

fprintf('Reset complete\n');
    
    %% Measurement Loop
    u2 = udp(remoteAddr, 'RemotePort', remotePort, 'LocalPort', localPort);
    try
        fopen(u2);
        connect=on;
        fprintf('Server connected and started\n');
    catch
        fclose(u2);
        connect=off;
        fprintf('Server not connected\n');
    end

  
res = inst.SendScpi('*IDN?');
assert(res.ErrCode == 0);
fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));

pause(0.01);
res = inst.SendScpi('*CLS'); % clear
assert(res.ErrCode == 0);

res = inst.SendScpi('*RST'); % reset
assert(res.ErrCode == 0);
%     sampleRateDAC = 9e9;
%     sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
%     res = inst.SendScpi(sampleRateDAC_str); % set sample clock
%     assert(res.ErrCode == 0);
fprintf('Reset complete\n');
% ---------------------------------------------------------------------
% ADC Config
% ---------------------------------------------------------------------
inst.SendScpi(':DIG:MODE DUAL');

inst.SendScpi(sprintf(':DIG:CHAN %d', adcChanInd)); 
inst.SendScpi(':DIG:ACQ:FREE');    
inst.SendScpi(sprintf(':DIG:FREQ %f', sampleRate));
inst.SendScpi(':DIG:DDC:CLKS AWG');
inst.SendScpi(':DIG:DDC:MODE COMP');
inst.SendScpi(':DIG:DDC:CFR2 75.38E6');
inst.SendScpi(':DIG:DDC:PHAS2 90');
inst.SendScpi(':DIG:CHAN:RANG LOW');
% Enable acquisition in the digitizer's channels  
inst.SendScpi(':DIG:CHAN:STATE ENAB');   

fprintf('ADC Configured\n');
fprintf('Clocks synced\n');
                
%% case 2: generating pulse sequence
amps = [1 1];
frequencies = [0 0];
lengths = [50.5e-6 50.5e-6];
fprintf("This is the length of the first pulse %d \n", lengths(1));
phases = [0 90];
mods = [0 0]; %0 = square, 1=gauss, 2=sech, 3=hermite 
spacings = [5e-6 43e-6];
markers = [1 1]; %always keep these on
markers2 = [0 0];
trigs = [0 1]; %acquire on every "pi" pulse

reps = [1 10000];
repeatSeq = [1]; % how many times to repeat the block of pulses

tof = -1000*(25.9874);

ch=1;
ch2 = 4;
initializeAWG(ch);
clearPulseDict();
clearBlockDict();

defPulse('init_pul', amps(1), mods(1), lengths(1), phases(1), spacings(1));
defPulse('theta1', amps(2), mods(2), lengths(2), phases(2), spacings(2));
defBlock('pulsed_SL', {'init_pul','theta1'}, reps(1:2), markers(1:2), trigs(1:2));
makeBlocks({'pulsed_SL'}, ch, repeatSeq);

% inst.SendScpi(sprintf(':INST:CHAN %d',ch2));
% inst.SendScpi(sprintf(':FREQ:RAST %d',2.5E9));
% %fprintf('Ch %s DAC clk freq %s\n', num2str(ch), num2str(sampleRateDAC)) 
% inst.SendScpi(':SOUR:VOLT MAX');
% inst.SendScpi(':INIT:CONT ON');
% 
% makeBlocks({'pulsed_SL'}, ch2, repeatSeq);
%generatePulseSeqIQ(ch, amps, frequencies, lengths, phases, mods, spacings, reps, markers, markers2, trigs);
%generatePulseSeqIQ(ch, amps, frequencies, lengths, phases, spacings, reps, markers, trigs, repeatSeq, indices);

setNCO_IQ(ch, 75.38e6+tof, 0);
inst.SendScpi(sprintf(':DIG:DDC:CFR2 %d', 75.38e6+tof));
fprintf('Calculate and set data structures...\n');
numberOfPulses_total = reps(2);
Tmax= 2;
tacq= 12;
numberOfPulses= floor(numberOfPulses_total/Tmax); %in 1 second %will be 1 for FID
loops=Tmax;

readLen = round2((tacq+2)*1e-6*2.7/(16*1e-9),96)-96;

offLen = 0;
rc = inst.SendScpi(sprintf(':DIG:ACQ:DEF %d, %d',numberOfPulses*loops, 2 * readLen));
assert(rc.ErrCode == 0);

inst.SendScpi(sprintf(':DIG:CHAN %d', adcChanInd))

rc = inst.SendScpi(':DIG:TRIG:SOUR EXT'); %digChan
assert(rc.ErrCode == 0);
%inst.SendScpi(':DIG:TRIG:TYPE GATE');
rc = inst.SendScpi(':DIG:TRIG:SLOP NEG');
assert(rc.ErrCode == 0)
rc = inst.SendScpi(':DIG:TRIG:LEV1 1.0');
assert(rc.ErrCode == 0)
rc = inst.SendScpi(sprintf(':DIG:TRIG:DEL:EXT %f', 12e-6)); % external trigger delay
assert(rc.ErrCode == 0)

fprintf('Instr setup complete and ready to aquire\n');

 rc = inst.SendScpi(':DIG:ACQ:FRAM:CAPT:ALL');   
 assert(rc.ErrCode == 0);
 rc = inst.SendScpi(':DIG:ACQ:ZERO:ALL');
 assert(rc.ErrCode == 0);

fprintf('Waiting... Listen for Shuttle\n');
rc = inst.SendScpi(':DIG:INIT OFF'); 
assert(rc.ErrCode == 0);
rc = inst.SendScpi(':DIG:INIT ON');
assert(rc.ErrCode == 0);

%% case 3: starting the pulse sequence
inst.SendScpi(sprintf(':DIG:CHAN 2'))
fprintf('Triggering pulse sequence\n');
rc = inst.SendScpi('*TRG');
                
    
    if cType == "LAN"
    inst.Disconnect();
    else
        admin.CloseInstrument(instId);    
        admin.Close();
    end     
clear inst;
clear;

%fclose();

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
    starttimeIdx = roundCheck(starttime, dt)+1;
    durationPts = roundCheck(pulseDuration, dt);

    if starttimeIdx > 0 && starttimeIdx + durationPts <= npoints
        rawSignal(starttimeIdx:starttimeIdx+durationPts) = sin(2*pi*freq*dt*(0:durationPts) + 2*pi*freq*starttime + phase);
        sine = ampScale(bits, rawSignal);
    else
        error('Pulse out of bounds');
    end
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
  dacSignal = (vertScaled + verticalScale);
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

function u = create_task_table(inst,task_list, awgSRate)
limit_cycles = 1.04e6;
task_list = arrange_tasktable(task_list,awgSRate,limit_cycles);
inst.SendScpi(sprintf(':TASK:COMP:LENG %d',length(task_list)));
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
        inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',segm)); %what is segment 3?? 
        inst.SendScpi(':TASK:COMP:TYPE SING');
        inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',cycles));
        
        if i==1 % define trigger
            inst.SendScpi(':TASK:COMP:ENAB TRG2');
        end
        
        if i==length(task_list)
            inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',1));
        else
            inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',i+1))
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

function task_list = build_tasktable(inst,pol_times,chirps,sampleRateDAC,varargin)
% function that create the task liste before creating the task table
% input: - pol_times : list that contains the different DNP time for each
% steps. 
%          - chirps: 2 dimension structure that contains the information
%          for each chirp, with positive and negative DNP respectively
%         
%          - sampleRateDAC: sample rate for the creation of the mw
%
% options : -  'first sign' could be '+' or '-' to stqrt with positive or
% negative DNP. positive DNP is the default value
%
%               - 'test' could be 'on' or 'off' depending on if the chirp
%               structure contains real chirps (on) or arbitrary wave form
%               for the test (off)
    
    test = 'off';
    chirp_num = 1;
    if isempty(varargin) == 0
        for i = length(varargin)/2
            if strcmp(varargin(i),'first sign') == 1
                if strcmp(varargin(i+1),'+') ==1
                    chirp_num = 1;
                elseif strcmp(varargin(i+1),'-')
                    chirp_num = 2;
                end
            elseif strcmp(varargin(i),'test') == 1
                test = varargin(i+1);
            end
        end
    end
    
    for i = 1:2
             % Define segment i 
        res = inst.SendScpi([':TRAC:DEF ',num2str(chirps{i}.segm),',',num2str(chirps{i}.segLen)]); % define memory location 1 of length dacSignal
        assert(res.ErrCode == 0);

        % select segment as the the programmable segment
        res = inst.SendScpi(sprintf(':TRAC:SEL %d',chirps{i}.segm));
        assert(res.ErrCode == 0); 

%         % Download the binary data to segment
%         res = inst.WriteBinaryData(':TRAC:DATA 0,', chirps{i}.dacSignal);
%         assert(res.ErrCode == 0);
        
        % Download the binary data to segment
        prefix = ':TRAC:DATA 0,';
        
        global bits
        if (bits==16)
            myWfm = uint16(chirps{i}.dacSignal);
            myWfm = typecast(myWfm, 'uint8');
        else
            myWfm = uint8(chirps{i}.dacSignal);
        end
        
        res = inst.WriteBinaryData(prefix, myWfm);
        
        if strcmp(test,'off') == 1
            srs_freq_str = [':SOUR:NCO:CFR1 ' sprintf('%0.2e',  sampleRateDAC - chirps{i}.srs_freq)]; %srs_freq
            res = inst.SendScpi(srs_freq_str);
            assert(res.ErrCode == 0);

           inst.SendScpi(':NCO:SIXD1 ON');

            rc = inst.SendScpi(':SOUR:MODE DUC');
            assert(rc.ErrCode == 0);
            
            try
                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);
                catch
                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);
            end
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
end

function initializeAWG(ch)
     global inst
     global interp
     global sampleRateInterp
     global sampleRateDAC
     
     sampleRateInterp = 5400e6;
     %sampleRateInterp = 9000e6
     %sampleRateInterp = 2017.5e6;
     %sampleRateInterp =  2*interp * sampleRateDAC;
     %sampleRateDAC = sampleRateDAC/(2*interp);
     sampleRateDAC = sampleRateInterp/(2*interp);
     inst.SendScpi(sprintf(':INST:CHAN %d',ch));
     %inst.SendScpi(sprintf(':FREQ:RAST %d',sampleRateDAC));
     inst.SendScpi(sprintf(':FREQ:RAST %d',2.5E9));
     %fprintf('Ch %s DAC clk freq %s\n', num2str(ch), num2str(sampleRateDAC)) 
     inst.SendScpi(':SOUR:VOLT MAX');
     inst.SendScpi(':INIT:CONT ON');
     res = inst.SendScpi(':TRAC:DEL:ALL');
     assert(res.ErrCode==0);
end



%function generatePulseSeqIQ(ch, amps, frequencies, lengths, phases, mods, spacings, reps, markers1, markers2, trigs)
function generatePulseSeqIQ(ch, amps, frequencies, lengths, phases, spacings, reps, markers1, trigs, repeatSeq, indices)
global inst
global sampleRateDAC
global interp
global sampleRateInterp
global blockDict
global pulseDict
    
    function downLoad_mrkr(ch, segMem, dacWave, mkrNum, state1, state2)
    fprintf('Downloading marker to channel %s, segment %s\n', num2str(ch), num2str(segMem))
    
    myMkr = uint8(state1 + 2*state2);
    
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi(sprintf(':TRAC:SEL %d',segMem));
    
    myMkr = myMkr(1:2:length(myMkr)) + 16 * myMkr(2:2:length(myMkr));

    res = inst.WriteBinaryData(':MARK:DATA 0,', myMkr);
    assert(res.ErrCode == 0);
    
    inst.SendScpi(sprintf(':MARK:SEL %d',1));
    inst.SendScpi(':MARK:VOLT:PTOP 0.5');
    %inst.SendScpi(':MARK:VOLT:LEV 0.0')
    inst.SendScpi(':MARK:VOLT:OFFS 0.25');
    inst.SendScpi(':MARK:STAT ON');
    
    inst.SendScpi(sprintf(':MARK:SEL %d',2));
    inst.SendScpi(':MARK:VOLT:PTOP 1.0');
    %inst.SendScpi(':MARK:VOLT:LEV 0.0')
    inst.SendScpi(':MARK:VOLT:OFFS 0.0');
    inst.SendScpi(':MARK:STAT ON');
    
    end

    function downLoadIQ(ch, segMem, dacWaveI, dacWaveQ, markerState1, markerState2, mkrNum)
        disp(sprintf('Downloading waveform to channel %s, segment %s', num2str(ch), num2str(segMem)))

        dacWaveIQ = [dacWaveI; dacWaveQ];
        dacWaveIQ = dacWaveIQ(:)';
        inst.SendScpi(sprintf(':INST:CHAN %d',ch));
        inst.SendScpi(':TRAC:FORM U16');
        inst.SendScpi(sprintf(':TRAC:DEF %d, %d',segMem, length(dacWaveIQ)));
        inst.SendScpi(sprintf(':TRAC:SEL %d',segMem));

%         res = inst.WriteBinaryData(':TRAC:DATA 0,', dacWaveIQ)
%         %assert(res.ErrCode==0);
        
        % Download the binary data to segment
        prefix = ':TRAC:DATA 0,';
        
        global bits
        if (bits==16)
            myWfm = uint16(dacWaveIQ);
            myWfm = typecast(myWfm, 'uint8');
        else
            myWfm = uint8(dacWaveIQ);
        end
        
        res = inst.WriteBinaryData(prefix, myWfm);
        
        downLoad_mrkr(ch, segMem, myWfm, mkrNum, markerState1, markerState2)
    end   

    function [mydcI, mydcQ] = makeDC(length)

    
    segLen = 64*round(length/64); %must be a multiple of 64
    
    max_dac = 2^16-1;
    half_dac = floor(max_dac/2);
    
    dacWave = zeros(1, segLen) + half_dac;
    
    mydcI = dacWave;
    mydcQ = dacWave;
    
    end 

    function [myWaveI, myWaveQ] = makeSqPulse(modFreq, pulseLen, amplitude, phase, mods)
        
        ampI = amplitude;
        ampQ = amplitude;

        segLen = 32*round(pulseLen/32); %must be a multiple of 32
        cycles = segLen * modFreq / sampleRateDAC;
        time = linspace(0, segLen-1, segLen);
        omega = 2 * pi * cycles;
    
        
        %disp('pulse modulation freq = ' + sampleRateDAC*cycles/segLen)
        if mods==1
            timeGauss = linspace(-segLen/2, segLen/2, segLen);
            sigma = segLen/6;
            modWave = exp(-0.5*(timeGauss/sigma).^2);

        elseif mods==2
            timeCosh = linspace(-segLen/2, segLen/2, segLen);
            tau = 2.355/1.76*segLen/6;
            modWave = cosh(timeCosh./tau).^-2;
        
        elseif mods==3
            timeHerm = linspace(-segLen/2, segLen/2, segLen);
            sigma = segLen/6;
            factor = 0.667;
            modWave = (1-factor*0.5*(timeHerm/sigma).^2).*exp(-0.5*(timeHerm/sigma).^2);
        
        else
            modWave = 1;
            
        end
        disp(sprintf('pulse segment length = %d points, actual time= %d', segLen, segLen/sampleRateDAC))
        max_dac = 2^16-1;
        half_dac = floor(max_dac/2);

        dacWave = ampI*cos(omega*time/segLen + pi*phase/180);
        dacWaveI = (dacWave.*modWave + 1)*half_dac;
        myWaveI = dacWaveI;
        
        dacWave = ampQ*sin(omega*time/segLen + pi*phase/180);
        dacWaveQ = (dacWave.*modWave + 1)*half_dac;
        myWaveQ = dacWaveQ;
    end 
    
    function setTask_Pulse(ch, numPulses, numSegs, reps, trigs, indices, repeatSeq)
    disp('setting task table')

    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi('TASK:ZERO:ALL');
    inst.SendScpi(sprintf(':TASK:COMP:LENG %d',numSegs)); % this should be more general?
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',1));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
    inst.SendScpi(':TASK:COMP:ENAB CPU');
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',1));
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',2));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    x=2;
    for y = 1:numBlocks
        lenBlock = length(indices{y});
        for z = 1:lenBlock
           inst.SendScpi(sprintf(':TASK:COMP:SEL %d',x));
           inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',indices{y}(z)));
           inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',reps{y}(z)));
           if (repeatSeq(y) > 1 && z==1) % if first task in a block
               inst.SendScpi('TASK:COMP:TYPE STAR');
               inst.SendScpi(sprintf(':TASK:COMP:SEQ %d',repeatSeq(y))); % number of loops for sequence   
           elseif (repeatSeq(y) > 1 && z~=lenBlock) % if intermediate task in a block
               inst.SendScpi('TASK:COMP:TYPE SEQ');
           elseif (repeatSeq(y) > 1 && z==lenBlock) % if last task in a block
               inst.SendScpi('TASK:COMP:TYPE END');
           else
               inst.SendScpi(':TASK:COMP:TYPE SING');
           end
           
           inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',x+1));
           
           if trigs{y}(z)==1
            inst.SendScpi('TASK:COMP:DTR ON');
           else
            inst.SendScpi('TASK:COMP:DTR OFF');
           end
           
           x=x+1;     
        end
    end
    
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',x));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',1));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',1));
    
    inst.SendScpi('TASK:COMP:WRITE');
    resp = inst.SendScpi('SOUR:FUNC:MODE TASK');
    assert(resp.ErrCode==0);
    end

    %%%% FUNCTION STARTS HERE %%%%
    numBlocks = length(blockDict);
    numPulses = 0;
    lengthsPts = {};
    spacingsPts = {};

    for u = 1:numBlocks
        numPulses = numPulses + length(amps{u});
        lengthsPts{end+1} = sampleRateDAC * lengths{u};
        spacingsPts{end+1} = sampleRateDAC * spacings{u};
    end
    numSegs = numPulses;
    disp('generating RF pulse sequence')
    global segMat
    segMat = cell(4, numSegs); % added a fourth row for Marker2 (trigs)

    %%%%%%% MAKE HOLDING SEGMENT %%%  
    DClen = 64;
    [holdI, holdQ] = makeDC(DClen);
    markHold = uint8(zeros(DClen, 1));
    
    x=1;
    for y = 1:numBlocks
        lenBlock = length(indices{y});
        for z = 1:lenBlock
            DClen = spacingsPts{y}(z);
            %%% make DC segment %%%
            [tempDCi, tempDCq] = makeDC(DClen);
            DClenreal = length(tempDCi);
            markDC = uint8(zeros(DClenreal, 1));
            markDC2 = uint8(zeros(DClenreal,1));
            %%% make Pulse %%% 
            [tempI, tempQ] = makeSqPulse(frequencies{y}(z), lengthsPts{y}(z),  amps{y}(z), phases{y}(z), 0);
            pulseLenReal = length(tempI);
            markIQ = uint8(zeros(pulseLenReal, 1) + markers1{y}(z));
            markIQ2 = uint8(zeros(pulseLenReal, 1) + trigs{y}(z));
            segMat{1,x} = [tempI tempDCi]; % first row is In-phase of pulse
            segMat{2,x} = [tempQ tempDCq]; % second row is Quadrature of pulse
            segMat{3,x} = [markIQ' markDC']; % third row is blanking signal (M1)
            segMat{4,x} = [markIQ2' markDC2']; % fourth row is ADC Ext trig signal (M2)
            x = x+1;
       end
    end
    %%% MAKE FINAL SEGMENT %%% 
    DClen = 64;
    [finalI, finalQ] = makeDC(DClen);

    fprintf('pulse sequence generated')
 
    downLoadIQ(ch, 1, holdI, holdQ, markHold, markHold, 1);
    x=1;
    for y = 1:numBlocks
        lenBlock = length(indices{y});
        for z = 1:lenBlock
           downLoadIQ(ch, indices{y}(z), segMat{1,x}, segMat{2,x}, segMat{3,x},segMat{4,x}, 1);
           x=x+1;
        end
    end
    downLoadIQ(ch, length(pulseDict)+2, finalI, finalQ, markHold, markHold, 1);
    setTask_Pulse(ch, numPulses, numSegs, reps, trigs, indices, repeatSeq);
 
    fprintf('pulse sequence written \n');
end

function setNCO_IQ(ch, cfr, phase)
    global sampleRateDAC
    global sampleRateInterp
    global inst
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi([':FREQ:RAST ' num2str(2.5E9)]);
    inst.SendScpi(':SOUR:INT X8');
    inst.SendScpi(':MODE DUC');
    inst.SendScpi(':IQM ONE');
    sampleRateDAC = sampleRateInterp;
    inst.SendScpi(sprintf(':FREQ:RAST %d', sampleRateDAC));
    inst.SendScpi('SOUR:NCO:SIXD1 ON');
    inst.SendScpi(sprintf(':SOUR:NCO:CFR1 %d',cfr));
    inst.SendScpi(sprintf(':SOUR:NCO:PHAS1 %d',phase));
    resp = inst.SendScpi(':OUTP ON');
    assert(resp.ErrCode==0);
end    


function outWfm = AlignTrig2(inWfm, rising_edge, pre_trigger, post_trigger)

    % inWmf always two channels. The first N acquisitions are for sync and the
    % last N acquisitions are for waveform to align.
    % rising_edge true = '+' edge, false = '-' edge.
    % pre_trigger, samples at the output before trigger.
    % post_trigger, samples at the ouput after trigger.
    % the real number of frames is half the number of frames in the input
    % array as it includes N frames for each channel (sync and waveform)
    % so total is 2 x N
    frame_length = size(inWfm, 2); % Just for information. Not used.
    num_of_frames = size(inWfm, 1);
    num_of_frames = num_of_frames / 2;
    % Threshold for edge detection is found as 50% of the Vlow to Vhigh of
    % the sync signal. Just the first frame is used to speed up the
    % process.
    threshold = (max(inWfm(1,:)) + min(inWfm(1,:))) / 2.0;

    % The aligned waveforms will be stored in this array
    % Actual length of the output frames is equal to the addition of
    % pre_trigger and post_trigger plus an additional sample for the
    % trigger itself.
    length_out = pre_trigger + post_trigger + 1;
    outWfm = zeros(num_of_frames,length_out);

    % The edge is found for all the sync signal frames
   
    for frame = 1:num_of_frames
        % If not found, the trigger position will be assigned to the
        % nominal trigger position.
        % The edge is searched in a range of samples around the nominal.
        % The width of the window must be higher than the maximum trigger
        % jitter to compensate (48 samples) plus any additional,
        % deterministic skew respect the nominal trigger position
        trig_sample = pre_trigger + 1;
        for sample = (pre_trigger + 1):frame_length
            % Threshold crossing is detected here
            if ((inWfm(frame, sample - 1) <= threshold && inWfm(frame, sample) > threshold) && rising_edge) ||...
                    ((inWfm(frame, sample - 1) >= threshold && inWfm(frame, sample) < threshold) && ~rising_edge)
                trig_sample = sample;
                break;
            end
        end
        % Time aligment is perfromed by selecting a portion of the
        % waveform to align around the effective trigger sample
        outWfm(frame, :) = inWfm(num_of_frames + frame,...
           (trig_sample - pre_trigger):(trig_sample + post_trigger));
    end
end

function defPulse(name, amp, mod, lengthTime, phase, waitTime)
    global pulseDict
    pulseName = name;
    index = length(pulseDict) + 2;
    pulseProps = [index, amp, mod, lengthTime, phase, waitTime, 0, 0, 0];
    pulseDict(pulseName) = pulseProps;
    fprintf('pulse %s added to dict at index %d\n', pulseName, index);
 end
       
function defBlock(blockName, pulseNames, reps, markers, trigs)
    global pulseDict
    global blockDict
    numPulses = length(pulseNames);
    blockProps = zeros(numPulses, 9);
    for x = 1:numPulses
        tempvec = pulseDict(pulseNames{x});
        tempvec(7) = reps(x);
        tempvec(8) = markers(x);
        tempvec(9) = trigs(x);
        pulseDict(pulseNames{x}) = tempvec;
        blockProps(x, :) = tempvec;
    end
    blockDict(blockName) = blockProps;
    fprintf('block %s created\n', blockName);
 end
 
function clearPulseDict()
    global pulseDict
    pulseDict = containers.Map;
end

function clearBlockDict()
    global blockDict
    blockDict = containers.Map;
end

function makeBlocks(blockNames, channel, repeatSeq)
    
    global blockDict
    
    indices = {};
    amps = {};
    frequencies = {};
    lengthsTime = {};
    phases = {};
    waitTimes = {};
    reps = {};
    markers = {};
    trigs = {};

    blockProps = 0;
    for x = 1:length(blockNames)
        name = blockNames{x}
        blockProps = blockDict(name);
        indices{end+1} = blockProps(:, 1);
        amps{end+1} = blockProps(:, 2);
        frequencies{end+1} = blockProps(:, 3);
        lengthsTime{end+1} = blockProps(:, 4);
        phases{end+1} = blockProps(:, 5);
        waitTimes{end+1} = blockProps(:, 6);
        reps{end+1} = blockProps(:, 7);
        markers{end+1} = blockProps(:, 8);
        trigs{end+1} = blockProps(:, 9);

    end
    disp('block made, passed to pulse gen')
    generatePulseSeqIQ(channel, amps, frequencies, lengthsTime, phases, waitTimes, reps, markers, trigs, repeatSeq, indices);
    
end
 
 