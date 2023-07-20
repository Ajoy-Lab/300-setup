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

% Define IP Address for Target Proteus device descriptor
% VISA "Socket-Based" TCP-IP Device. Socket# = 5025
ipAddr = '127.0.0.1'; %'127.0.0.1'= Local Host; % your IP here
pxiSlot = 0;

% Instrument setup
cType =         "DLL";  %"LAN" = VISA or "DLL" = PXI

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

channel = 1;
% AWG Settings
samplingRate = 2500E6;
interpol     = 1;
dacMode = 16;
wfmVolt = 0.25;
wfmOff = 0.0;

% Marker Settings
marker_pos = 150E-9;
marker_width = 150E-9;
mkrVolt = 1.0;
mkrOff = 0.5;

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

% Segment
numOfSegments = 4;

fprintf(1, 'Calculating WAVEFORMS\n');

minCycles = 1;
period = 1E-7;

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

%Select Channel
inst.SendScpi(sprintf(':INST:CHAN %d', channel));
% DAC Mode set to 'DIRECT" (Default)
inst.SendScpi(':SOUR:MODE DIRECT');

% Set sampling rate for AWG to maximum.
if interpol > 1
    inst.SendScpi(':FREQ:RAST 2.5E9');
    inst.SendScpi([':INT X ' num2str(interpol)]);
end
inst.SendScpi([':FREQ:RAST ' num2str(samplingRate)]);

for segment = 1:numOfSegments
    % Calculate basic square wave
    if mod(segment, 4) == 1
        myWfm = getSquareWfm(   samplingRate / interpol,... 
                                minCycles,...
                                period,...
                                granul);
    end

    mkrDiv = 2;
    if dacMode == 8
        mkrDiv = 8;
    end

    % Marker processing

    mkr_sampling_rate = samplingRate / mkrDiv;
    mkr_length = length(myWfm) / mkrDiv;
    
    start_sample_marker = round(marker_pos * mkr_sampling_rate);
    stop_sample_marker = round((marker_pos + marker_width) * mkr_sampling_rate);
    
    total_length = numOfSegments * mkr_length;

    if start_sample_marker > total_length
        start_sample_marker = 0; % No marker generation
    end

    if stop_sample_marker > total_length
        stop_sample_marker = total_length; % No marker generation
    end

    myMkr1 = uint8(zeros(1, total_length));
    if start_sample_marker > 0
        myMkr1(start_sample_marker:stop_sample_marker) = uint8(1);  
        myMkr2 = uint8(zeros(1, mkr_length));         
        % Select portion
        myMkr1 = myMkr1(((segment - 1) * mkr_length + 1):(segment * mkr_length));        
    else
        myMkr1 = uint8(zeros(1, mkr_length));
        myMkr2 = uint8(zeros(1, mkr_length));
    end
    
    myMkr = FormatMkr2(dacMode, myMkr1, myMkr2);
    
    % Segment # processing
    % All Proteus models except the P908X share the same waveform memory
    % bank among channel N+1 and N+2, N=0..NumOfChannels/2. This means that
    % the same segment number cannot be used for this pair of channels. In
    % this case the designated segment is used for the odd numbered
    % channels and the next segment is assigned to the even numbered channel
    % of the same pair. All segments can be deleted just once for each pair
    % of channels.
    if segment == 1     
        % All segments deleted for current waveform memory bank
        inst.SendScpi(':TRAC:DEL:ALL');
    end    
    
    % Waveform Downloading
    % *******************
    
    fprintf(1, 'DOWNLOADING WAVEFORM FOR SEGMENT%d\n', segment);

    SendWfmToProteus(       inst,...
                            samplingRate,...
                            channel,...
                            segment,...
                            myWfm,...
                            dacMode,...
                            false);

    result = SendMkrToProteus(inst, myMkr);

    fprintf(1, 'WAVEFORM DOWNLOADED!\n');
    % Select segment for generation
    fprintf(1, 'SETTING AWG OUTPUT\n');

    inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segment));    
    
    % The new waveform is calculated for the next channel
    if mod(segment - 1, 4) < 3
        % Integration
        myWfm = cumsum(myWfm);
        % DC removal
        myWfm = myWfm - mean(myWfm);
        % Normalization to the -1.0/+1.0 range
        myWfm = myWfm / max(abs(myWfm));
    end
end

% TASK LIST GO HERE

inst.SendScpi(sprintf(':TASK:COMP:LENG %d', numOfSegments));

for k = 1:numOfSegments
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d', k));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    inst.SendScpi(':TASK:COMP:DEST NEXT');
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', k));
    nexTask = k + 1;
    if k == numOfSegments
        nexTask = 1;
    end
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d', nexTask));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', 1));
end

inst.SendScpi(':TASK:COMP:WRIT');
fprintf(1, 'SEQUENCE CREATED!\n');
% Select Task Mode for generation    
% Start in task #1 (#1 is the default)    
inst.SendScpi(':FUNC:MODE TASK');


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

% It is recommended to disconnect from instrument at the end
if cType == "LAN"
    inst.Disconnect();
else
    admin.CloseInstrument(inst.InstrId);    
    admin.Close();
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


