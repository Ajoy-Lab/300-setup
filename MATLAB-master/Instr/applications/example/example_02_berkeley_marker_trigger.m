% EXAMPLE FOR DIRECT MODE
%===================================================
% This example calculates up to 4 different signals and download them into
% each available channel in the target Proteus device.
% 
% The basic waveform is an square waveform using the full DAC range and it
% is downloaded to channel #1. For each channel, the waveform is calculated
% by integration of the previous waveform in a similar way to some analog
% signal generators, where the triangular wave is obtained by integration
% of an square wave, and the sinusoidal waveform is obtained by integration
% of the triangular wave. Channel #4, when available, will generate a
% "cosine" wave obatined by integration of the sinewave assigned to channel
% #3.
% Markers for each channel are also calculated and downloaded. Marker 1 is
% a sync pulse. Its duration (in states) is equal to the channel number.
% Marker 2 is just geenrating a random stream of bits.

clear;
close all;
clear variables;
clear global;
clc;

fprintf(1, 'INITIALIZING SETTINGS\n');

pid = feature('getpid');
fprintf(1,'\nProcess ID %d\n',pid);

% BASIC EXAMPLE FOR CONNECTION TO PROTEUS USING VISA OR PXI
%==========================================================
% VISA Communications from MATLAB requires the Instrument Control Toolbox

clear;
close all;
clear variables;
clear global;
clc;

% Define IP Address for Target Proteus device descriptor
% VISA "Socket-Based" TCP-IP Device. Socket# = 5025
ipAddr = '127.0.0.1'; %'127.0.0.1'= Local Host; % your IP here
pxiSlot = 0;

% Instrument setup
cType =         "LAN";  %"LAN" = VISA or "DLL" = PXI

if cType == "LAN"
    connPar =       ipAddr; 
else
    connPar =       pxiSlot; % Your slot # here, o for manual selection
end

paranoia_level = 2; % 0, 1 or 2
% Open Session and load libraries
[inst, admin, model, slotNumber] = ConnecToProteus(cType, connPar, paranoia_level);

% Report model
fprintf('Connected to: %s, slot: %d\n', model(1), slotNumber(1));

% Reset AWG
inst.SendScpi('*CLS;*RST');

% Get options using the standard IEEE-488.2 Command
optstr = getOptions(inst);

samplingRate = 2500E6;
interpol     = 1;
dacMode = 16;

if (samplingRate / interpol) > 2.5E9
    dacMode = 8;
    interpol = 1;
end

if (samplingRate / interpol) < 250E6
    interpol = 1;
end

% Get granularity
granul = getGranularity(model(1), optstr, dacMode);
fprintf('\nGranularity = %d samples\n', granul);

% Get Active Channels and Segment #
[chanList, segmList] = GetChannels(model(1), samplingRate / interpol);
numOfChannels = length(chanList);

fprintf(1, 'Calculating WAVEFORMS\n');

minCycles = 1;
period = 1E-7;

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Set sampling rate for AWG to maximum.
if interpol > 1
    inst.SendScpi(':FREQ:RAST 2.5E9');
    inst.SendScpi([':INT X ' num2str(interpol)]);
end
inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);

wfmVolt = 0.25;
wfmOff = 0.0;

mkrVolt = 0.1;
mkrOff = 0.0;

for channel = 1:numOfChannels
    % Calculate basic square wave
    if mod(channel, 4) == 1
        myWfm = getSquareWfm(   samplingRate / interpol,... 
                                minCycles,...
                                period,...
                                granul);
    end

    mkrDiv = 2;
    if dacMode == 8
        mkrDiv = 8;
    end

    myMkr1 = uint8(zeros(1, length(myWfm) / mkrDiv));
    myMkr1(1:40) = uint8(1); % Sync marker duration depends on the channel
    myMkr2 = rand(1, length(myMkr1)); % Random data is different for each channel
    myMkr2 = uint8(myMkr2 > 0.5);
    myMkr = FormatMkr2(dacMode, myMkr1, myMkr2);
%     myMkr2 = uint8(myMkr2 > 0.5);
%     myMkr = myMkr + 2 * myMkr2;
% 
%     if dacMode == 16
%         myMkr = myMkr(1:2:length(myMkr)) + 16 * myMkr(2:2:length(myMkr));
%     end

    %Select Channel
    inst.SendScpi(sprintf(':INST:CHAN %d', chanList(channel)));
    % DAC Mode set to 'DIRECT" (Default)
    inst.SendScpi(':SOUR:MODE DIRECT');
    
    % Segment # processing
    % All Proteus models except the P908X share the same waveform memory
    % bank among channel N+1 and N+2, N=0..NumOfChannels/2. This means that
    % the same segment number cannot be used for this pair of channels. In
    % this case the designated segment is used for the odd numbered
    % channels and the next segment is assigned to the even numbered channel
    % of the same pair. All segments can be deleted just once for each pair
    % of channels.
    if segmList(channel) == 1     
        % All segments deleted for current waveform memory bank
        inst.SendScpi(':TRAC:DEL:ALL');
    end    
    
    % Waveform Downloading
    % *******************
    
    fprintf(1, 'DOWNLOADING WAVEFORM FOR CH%d\n', chanList(channel));

    SendWfmToProteus(       inst,...
                            samplingRate,...
                            chanList(channel),...
                            segmList(channel),...
                            myWfm,...
                            dacMode,...
                            false);

    result = SendMkrToProteus(inst, myMkr);

    fprintf(1, 'WAVEFORM DOWNLOADED!\n');
    % Select segment for generation
    fprintf(1, 'SETTING AWG OUTPUT\n');

    inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segmList(channel)));

    % Output voltage and offset
    inst.SendScpi([':SOUR:VOLT ' num2str(wfmVolt)]);  
    inst.SendScpi([':SOUR:VOLT:OFFS ' num2str(wfmOff)]); 
    % Activate outpurt and start generation
    inst.SendScpi(':OUTP ON');   

    inst.SendScpi(':MARK:SEL 1');         %Marker1
    inst.SendScpi([':MARK:VOLT:PTOP ' num2str(mkrVolt)]);      %Vpp
    inst.SendScpi([':MARK:VOLT:OFFS ' num2str(mkrOff)]);      %DC Offset   
    inst.SendScpi(':MARK ON');   
    inst.SendScpi(':MARK:SEL 2');         %Marker1
    inst.SendScpi([':MARK:VOLT:PTOP ' num2str(mkrVolt)]);      %Vpp
    inst.SendScpi([':MARK:VOLT:OFFS ' num2str(mkrOff)]);      %DC Offset 
    inst.SendScpi(':MARK ON');
    
    % The new waveform is calculated for the next channel
    if mod(channel - 1, 4) < 3
        % Integration
        myWfm = cumsum(myWfm);
        % DC removal
        myWfm = myWfm - mean(myWfm);
        % Normalization to the -1.0/+1.0 range
        myWfm = myWfm / max(abs(myWfm));
    end
end

% Get all frames from digitizer
fprintf('Acquired Waveforms Upload to Computer Starts\n');

% ***************************************************************
% Relevant Variables for ADC in Dual Mode, no DDC, self trigger
dualMode = true;
digChannel = [2, 1];
adcRange = [1, 2];
use_ddc = false;
use_dtr_trigger = false;
trgLevel = 0.0;
adc_sampling_rate = 2.7E9;
frame_length = 960;
num_of_frames = 100;
pre_trigger = 480;

% Non-Relevant Variables
awg_channel = 1;
dtrDelay = 0.0;
carrier_freq = 0.0E9;


acqWfm = GetDigitizerData(  inst,...
                            cType,...
                            dualMode,...
                            digChannel,...
                            adcRange,...
                            use_ddc,...
                            use_dtr_trigger,...
                            awg_channel,...
                            dtrDelay,...
                            trgLevel,...
                            adc_sampling_rate,...
                            carrier_freq,...                              
                            frame_length,...
                            num_of_frames,...
                            pre_trigger);

outWfm = AlignTrig(acqWfm, pre_trigger, frame_length - 200);

tiledlayout(1,3);
ax1 = nexttile;

for frame = 1:num_of_frames
    plot(ax1, acqWfm(frame,:));
    if frame == 1
        hold;
    end
end

ax2 = nexttile;
for frame = (num_of_frames + 1):(2 * num_of_frames)
    plot(ax2, acqWfm(frame,:));
    if frame == (num_of_frames + 1)
        hold;
    end
end

ax3 = nexttile;

for frame = 1:num_of_frames
    plot(ax3, outWfm(frame,:));
    if frame == 1
        hold;
    end
end


fprintf('Acquired Waveforms Upload to Computer Ends\n'); 

% It is recommended to disconnect from instrument at the end
if cType == "LAN"
    inst.Disconnect();
else
    admin.CloseInstrument(inst.InstrId);    
    admin.Close();
end  

function outWfm = AlignTrig(inWfm, trig_pos, length_out)

    % inWmf always two channels
    % trig_pos in samples
    % length_out, samples at the output
    % trigger signal is alays the first one

    frame_length = size(inWfm, 2);
    num_of_frames = size(inWfm, 1);
    num_of_frames = num_of_frames / 2;

    threshold = (max(inWfm(1,:)) + min(inWfm(1,:))) / 2.0; 

    
    outWfm = zeros(num_of_frames,length_out + 1);
    for frame = 1:num_of_frames
        trig_sample = trig_pos;
        for sample = (trig_pos - 60):(trig_pos + 60)
            if inWfm(frame, sample - 1) <= threshold && inWfm(frame, sample) > threshold
                trig_sample = sample;
                break;
            end
        end

        outWfm(frame, :) = inWfm(num_of_frames + frame,...
            (trig_sample - length_out/ 2):(trig_sample + length_out/ 2));
    end
end

function wfmData = GetDigitizerData(inst,...
                                    cType,...
                                    dualMode,...
                                    adChan,...
                                    range,...
                                    useDdc,...
                                    useDtrTrig,...                                    
                                    dtrChan,...
                                    dtrDelay,...
                                    trgLevel,...
                                    samplingRateDig,...
                                    cFreq,...
                                    frameLen,...
                                    numberOfFrames,...
                                    preTrig)
    % ----------
    % ADC Config
    % It supports up to two channels with different setups
    % ----------

    % Set the ADC mode and set the channel mapping
    if dualMode
        inst.SendScpi(':DIG:MODE DUAL');
        if length(adChan) > 1
            adChan = adChan(1:2);
            if adChan(1) == adChan(2)
                adChan = adChan(1);
            end
        end        
    else
        inst.SendScpi(':DIG:MODE SINGLE');
        % Only Channel 1 in Single mode
        adChan = 1;        
    end  

    adChan(adChan > 2) = 2;
    adChan(adChan < 1) = 1;

    numOfChannels = length(adChan);

    % Free Acquistion Memory and Set sampling rate
    inst.SendScpi(':DIG:ACQ:FREE');    
    inst.SendScpi(sprintf(':DIG:FREQ %g', samplingRateDig));
    % DDC activation
    if useDdc        
        inst.SendScpi(':DIG:DDC:MODE COMP');   
        for i = 1:numOfChannels
            % NCO frequency mapped to the first and second NZ
            ddcFreq = GetNcoFreq(cFreq(mod(i, length(cFreq)) + 1), samplingRateDig, false);
            inst.SendScpi([sprintf(':DIG:DDC:CFR%d ', adChan(i)) num2str(abs(ddcFreq))]);
        end
    end
    % Calculate actual frame length depending on the DDC mode
    actualFrameLen = frameLen;
    if useDdc
        actualFrameLen = 2 * actualFrameLen;
    end
    
    % ADC Range 
    % 1:LOW, 2:MED, 3:HIGH
    range(range > 2) = 3;
    range(range < 2) = 1;
    
    for chan = 1:numOfChannels
        % Select digitizer channel:
        inst.SendScpi(sprintf(':DIG:CHAN %d', adChan(chan)));
        % Set the voltage-range of the selected channel
        switch range(mod(chan - 1, length(range)) + 1)
            case 1
                inst.SendScpi(':DIG:CHAN:RANG LOW');
            case 2
                inst.SendScpi(':DIG:CHAN:RANG MED');
            case 3
                inst.SendScpi(':DIG:CHAN:RANG HIGH');
        end
        %Enable acquisition in the selected channel
        inst.SendScpi(':DIG:CHAN:STATE ENAB');

        % Setup frames layout. Common to both ADC channels.    
        inst.SendScpi(sprintf(':DIG:ACQ:DEF %d, %d',...
            numberOfFrames, actualFrameLen));
        
        % Set channel 1 of the digitizer as its trigger source
        % If DTR trigger, it is directed to the designated AWG channel
        if useDtrTrig(mod(chan - 1, length(useDtrTrig)) + 1) 
            % DTR trigger must be assigned after selecting the target AWG
            % channel as it is a property of the AWG channel and the ADC
            % channel
            inst.SendScpi(sprintf(':INST:CHAN %d',...
                dtrChan(mod(chan - 1, length(dtrChan)) + 1)));
            inst.SendScpi(sprintf(':DIG:TRIG:SOURCE TASK%d',...
                dtrChan(mod(chan - 1, length(dtrChan)) + 1)));
            % Set DTR trigger Dealy
            inst.SendScpi(sprintf(':DIG:TRIG:AWG:TDEL %f',...
                dtrDelay));
        else
            % Level trigger set tup
            inst.SendScpi(sprintf(':DIG:TRIG:SOUR CH%d', adChan(1))); %i
            if chan == adChan(1)
                inst.SendScpi(sprintf(':DIG:TRIG:SELF %f',...
                    trgLevel(mod(chan - 1, length(trgLevel)) + 1))); %0.025
            end
        end
        % Pretrigger for DDC must be set to double as acquisions are made
        % by IQ pair of samples
        if useDdc
            actualPreTrig = 2 * preTrig;
        else
            actualPreTrig = preTrig;
        end
    
        inst.SendScpi(sprintf(':DIG:PRET %d', actualPreTrig));
        
        % Select which frames are filled with captured data 
        %(all frames in this example)
        inst.SendScpi(':DIG:ACQ:FRAM:CAPT:ALL');
    
        % Delete all wfm memory
        inst.SendScpi(':DIG:ACQ:ZERO:ALL');
        % Get ADC wfm format. For informative purposes
        resp = inst.SendScpi(':DIG:DATA:FORM?');
        resp = strtrim(netStrToStr(resp.RespStr));         
    end
    
    % Stop the digitizer
    inst.SendScpi(':DIG:INIT OFF');
    % And start for a new acquisition
    wfmData = zeros(numOfChannels, frameLen);
    inst.SendScpi(':DIG:INIT ON'); 
    % Read Acquired Wfm Data
    for i=1:numOfChannels             
        % Select channel
        inst.SendScpi(sprintf(':DIG:CHAN %d', adChan(i)));
        % Get acquisition status CSV string from Proteus for selected
        % channel
        for n = 1:250
            resp = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
            resp = strtrim(netStrToStr(resp.RespStr));
            resp = strtrim(resp);
            items = split(resp, ',');
            items = str2double(items);
            % If item 2 in the CSV string is '1', then all frames have been
            % captured
            if length(items) >= 3 && items(2) == 1
                break
            end
            % This is just to give some information when trigger times out
            if mod(n, 10) == 0                
                fprintf('%d. %s Time:\n', fix(n / 10), resp);                                
            end
            pause(0.1);
        end          

        % Define what we want to read 
        % (frames data, frame-header, or both).
        % In this example we read the frames-data
        inst.SendScpi(':DIG:DATA:TYPE FRAM');
        inst.SendScpi(':DIG:DATA:SEL ALL');

        % Read binary block
        % Get the size in bytes of the acquisition
        resp = inst.SendScpi(':DIG:DATA:SIZE?');
        resp = strtrim(netStrToStr(resp.RespStr));
        num_bytes = str2double(resp);
        % upload time will be shown so transfer rate can be compared
        fprintf('ADC Upload Time for %d bytes:\n', num_bytes);
        
        if cType == "LAN"
            if useDdc 
                tic;
                % DDC data is formatted as 15-bit in a 32-bit unsigned
                % integer
                samples = inst.ReadBinaryData(':DIG:DATA:READ?', 'uint32');              
                toc;
                samples = int32(samples) - 16384; % Set zero level
                % Convert to complex I + jQ samples
                samples = InterleavedToComplex(samples);
                % Invert spectrum if ddcFreq < 0.0
                if ddcFreq < 0.0
                    samples = conj(samples);
                end
            else
                tic;
                % Direct ADC data is formated as 12-bit samples in 16-bit
                % unsigned integers
                samples = inst.ReadBinaryData(':DIG:DATA:READ?', 'uint16');     
                toc;
                samples = int16(samples) - 2048; % Set zero level
            end
        else
            % For the PXI library, downloads can only handlw 8 or 16-bit
            % unsigned.
            if useDdc
                % For DDC, because read format is UINT16 we divide byte
                % number by 2
                wavlen = floor(num_bytes / 2);        
                % allocate NET array
                netArray = NET.createArray('System.UInt16', wavlen);
                % read the captured frame
                tic;
                res = inst.ReadMultipleAdcFrames(i - 1, 1, numberOfFrames, netArray);
                toc;
             
                assert(res == 0);
                % Each 32 sample is now 2 contiguous 16-bit samples
                samples = uint16(netArray);
                % As the first 16-bit samples in the pair is "all zeros"
                % they can be discarded by taking one very two bytes
                samples = samples(1:2:length(samples));
                
                % cast to matlab vector
                samples = int16(samples) - 16384; % Set zero level
                % Convert to complex I + jQ samples
                samples = InterleavedToComplex(samples);
                % Invert spectrum if ddcFreq < 0.0
                if ddcFreq < 0.0
                    samples = conj(samples);
                end                
                % deallocate the NET array
                delete(netArray);
            else
                wavlen = floor(num_bytes / 2);
        
                % allocate NET array
                netArray = NET.createArray('System.UInt16', wavlen);
                % read the captured frame
                tic
                res = inst.ReadMultipleAdcFrames(0, 1, numberOfFrames, netArray);
                toc
                assert(res == 0);
                
                samples = uint16(netArray);
                % cast to matlab vector
                samples = int16(samples) - 2048;
                                
                % deallocate the NET array
                delete(netArray);
            end            
        end        
        % Ouput data is formatted as a two dimensions array with A x F
        % rows (A = number of acquisitions, F = number of Frames) and
        % FrameLem columns
        for j=1:numberOfFrames
            wfmData((i - 1) * numberOfFrames + j,:) = samples(((j-1) * frameLen + 1):(j * frameLen));
        end        
    end
    % Sttop digitizer after all acquisitions and frames for all the
    % channels have been captured
    inst.SendScpi(':DIG:INIT OFF');
end

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

function result = SendWfmToProteus(     inst,...
                                        samplingRate,...
                                        channel,...
                                        segment,...
                                        myWfm,...
                                        dacRes,...
                                        initialize)

    if dacRes == 16  
            inst.SendScpi(':TRAC:FORM U16');
    else
            inst.SendScpi(':TRAC:FORM U8');
    end

    %Select Channel
    if initialize
        inst.SendScpi(':TRAC:DEL:ALL');
        inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);        
    end
    
    inst.SendScpi(sprintf(':INST:CHAN %d', channel));    
    inst.SendScpi(sprintf(':TRAC:DEF %d, %d', segment, length(myWfm)));        
    % select segmen as the the programmable segment
    inst.SendScpi(sprintf(':TRAC:SEL %d', segment));

    % format Wfm
    myWfm = myQuantization(myWfm, dacRes, 0);
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';

    if (dacRes==16)
        myWfm = uint16(myWfm);
        myWfm = typecast(myWfm, 'uint8');
    else
        myWfm = uint8(myWfm);
    end
    
    res = inst.WriteBinaryData(prefix, myWfm);    
    assert(res.ErrCode == 0);

    if initialize
        inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segment))
        % Output voltage set to MAX
        inst.SendScpi(':SOUR:VOLT MAX');   
        % Activate outpurt and start generation
        inst.SendScpi(':OUTP ON');        
    end
    
    result = length(myWfm);
end 

function result = SendMkrToProteus( inst,myMkr)
    % Download the binary data to segment   
    prefix = ':MARK:DATA 0,';
    inst.WriteBinaryData(prefix, myMkr); 
    %instHandle.SendBinaryData(prefix, myMkr, 'uint8');
    result = length(myMkr);
end

function mkrData = FormatMkr2(dac_Mode, mkr1, mkr2)
    % Mkr1 goes to bit 0 and Mkr2 goes to bit 1 in a 4-bit Nibble
    mkrData = mkr1 + 2 * mkr2;
    % For DAC Mode 8, just one Nibble per Byte is sent
    % For DAC Mode 16, two consecutive nibbles are multiplexed in one byte
    if dac_Mode == 16
        mkrData = mkrData(1:2:length(mkrData)) + ...
            16 * mkrData(2:2:length(mkrData));
    end
end

function mkrData = FormatMkr4(dac_Mode, mkr1, mkr2, mkr3, mkr4)
    % Mkr1 goes to bit 0 and Mkr2 goes to bit 1 in a 4-bit Nibble
    mkrData = mkr1 + 2 * mkr2 + 4 * mkr3 + 8 * mkr4;
    % For DAC Mode 8, just one Nibble per Byte is sent
    % For DAC Mode 16, two consecutive nibbles are multiplexed in one byte
    if dac_Mode == 16
        mkrData = mkrData(1:2:length(mkrData)) + ...
            16 * mkrData(2:2:length(mkrData));
    end
end

function [  inst,...
            admin,...
            modelName,...
            sId] = ConnecToProteus( cType, ...
                                    connStr, ...
                                    paranoia_level)

% Connection to target Proteus
% cType specifies API. "LAN" for VISA, "DLL" for PXI
% connStr is the slot # as an integer(0 for manual selection) or IP adress
% as an string
% Paranoia Level add additional checks for each transfer. 0 = no checks.
% 1 = send OPC?, 2 = send SYST:ERROR?

% It returns
% inst: handler for the selected instrument
% admin: administrative handler
% modelName: string with model name for selected instrument (i.e. "P9484")
% sId: slot number for selected instrument   
    
    pid = feature('getpid');
    fprintf(1,'\nProcess ID %d\n',pid);
    
    dll_path = 'C:\\Windows\\System32\\TEPAdmin.dll';  
    admin = 0;

    sId = 0;
    
    if cType == "LAN"
        try            
            connStr = strcat('TCPIP::',connStr,'::5025::SOCKET');
            inst = TEProteusInst(connStr, paranoia_level);
            
            res = inst.Connect();
            assert (res == true);
            [modelName, sId] = identifyModel(inst);
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
            
            % If there are multiple slots, let the user selecidentifyModelt one ..
            sId = slotIds(1);
            if numSlots > 1
                fprintf('\n%d slots were found\n', numSlots);
                for n = 1:numSlots
                    sId = slotIds(n);
                    slotInfo = admin.GetSlotInfo(sId);
                    if ~slotInfo.IsSlotInUse
                        modelName = slotInfo.ModelName;
                        if slotInfo.IsDummySlot && connStr == 0
                            fprintf(' * Slot Number:%d Model %s [Dummy Slot].\n', sId, modelName);
                        elseif connStr == 0
                            fprintf(' * Slot Number:%d Model %s.\n', sId, modelName);
                        end
                    end
                end
                pause(0.1);
                if connStr == 0
                    choice = input('Enter SlotId ');
                    fprintf('\n');
                else
                    choice = connStr;
                end                
                sId = uint32(choice);
                slotInfo = admin.GetSlotInfo(sId);
                modelName = slotInfo.ModelName;
                modelName = strtrim(netStrToStr(modelName));
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
end

function [model, slot] = identifyModel(inst)
    idnStr = inst.SendScpi('*IDN?');
    idnStr = strtrim(netStrToStr(idnStr.RespStr));
    idnStr = split(idnStr, [",", "-", " "]);

    model = string.empty;
    slot = uint8.empty;

    for k = 1:length(idnStr)
        if contains(idnStr(k), "P")
           model = [model, idnStr(k)];
        end
        if contains(idnStr(k), "slot")
           slot = [slot, uint8(str2num(idnStr(k + 1)))];
        end
    end    
end

function options = getOptions(inst)
    optStr = inst.SendScpi('*OPT?');  
    optStr = strtrim(netStrToStr(optStr.RespStr));
    options = split(optStr, ',');  
end

function granularity = getGranularity(model, options, dacMode)

    flagLowGranularity = false;

    for i = 1:length(options)
        if contains(options(i), 'G1') || contains(options(i), 'G2')
            flagLowGranularity = true;
        end        
    end   

    flagLowGranularity = false; % TEMPORARY

    granularity = 32;

    if contains(model, 'P258')
        granularity = 32;    
        if flagLowGranularity
            granularity = 16;
        end
    elseif contains(model, 'P128')
        granularity = 32;    
        if flagLowGranularity
            granularity = 16;
        end
    elseif contains(model, 'P948')
        if dacMode == 16
            granularity = 32;    
            if flagLowGranularity
                granularity = 16;
            end
        else
            granularity = 64;    
            if flagLowGranularity
                granularity = 32;
            end
        end
    elseif contains(model, 'P908')
        granularity = 64;    
        if flagLowGranularity
            granularity = 32;
        end
    end
end

function [chanList, segmList] = GetChannels(model, sampleRate)

    if contains(model, 'P9484') || contains(model, 'P2584') || contains(model, 'P1284')
        if sampleRate <= 2.5E9
            chanList = [1 2 3 4];
            segmList = [1 2 1 2];
        else
            chanList = [1 3];
            segmList = [1 1];
        end

    elseif contains(model, 'P9482') || contains(model, 'P2582') || contains(model, 'P1282')
        if sampleRate <= 2.5E9
            chanList = [1 2];
            segmList = [1 2];
        else
            chanList = [1];
            segmList = [1];
        end

    elseif contains(model, 'P9488') || contains(model, 'P2588') || contains(model, 'P1288')
        if sampleRate <= 2.5E9
            chanList = [1 2 3 4 5 6 7 8];
            segmList = [1 2 1 2 1 2 1 2];
        else
            chanList = [1 3 5 7];
            segmList = [1 1 1 1];
        end

    elseif contains(model, 'P94812') || contains(model, 'P25812') || contains(model, 'P12812')
        if sampleRate <= 2.5E9
            chanList = [1 2 3 4 5 6 7 8 9 10 11 12];
            segmList = [1 2 1 2 1 2 1 2 1 2 1 2];
        else
            chanList = [1 3 5 7 9 11];
            segmList = [1 1 1 1 1 1];
        end

    elseif contains(model, 'P9082')        
        chanList = [1 2];
        segmList = [1 1];

    elseif contains(model, 'P9084')        
        chanList = [1 2 3 4];
        segmList = [1 1 1 1];

    elseif contains(model, 'P9086')        
        chanList = [1 2 3 4 5 6];
        segmList = [1 1 1 1 1 1];
    end
end

function retval = myQuantization (myArray, dacRes, minLevel)  
 
    maxLevel = 2 ^ dacRes - 1;  
    numOfLevels = maxLevel - minLevel + 1;
    
    retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
    retval = retval + minLevel;
    
    retval(retval > maxLevel) = maxLevel;
    retval(retval < minLevel) = minLevel;

end

function [str] = netStrToStr(netStr)
    try
        str = convertCharsToStrings(char(netStr));
    catch        
        str = '';
    end
end


