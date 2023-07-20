
% EXAMPLE FOR DIRECT MODE
%===================================================
% % This example calculates up to 4 different signals and download them into
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

clear;
close;
clc;

pfunc = ProteusFunctions;

pid = feature('getpid');
fprintf(1,'\nProcess ID %d\n',pid);

dll_path = 'C:\\Windows\\System32\\TEPAdmin.dll';

% Communication Parameters
connStr =           '127.0.0.1'; % your IP here
paranoia_level =    2; % 0, 1 or 2

cType = "DLL";  %"LAN" or "DLL"

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
            choice = input('Enter SlotId ');
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

try
    %Get instrument ID
    res = inst.SendScpi('*IDN?');
    fprintf('Connected to: %s\n', pfunc.netStrToStr(res.RespStr));
    
    %get number of channels
    res=inst.SendScpi(':INST:CHAN? MAX');
    num_channels=str2double(pfunc.netStrToStr(res.RespStr));
    
    % Get Model name
    res=inst.SendScpi(':SYST:INF:MODel?');
    model_name = pfunc.netStrToStr(res.RespStr);
    
    % Setup model dependant parameters 
    if strncmp(model_name,'P908',4)
        dac_mode = 8;
        num_ddr = num_channels;
        markers_per_chan = 4;
        max_dac = 255;
        data_type = 'uint8';
    elseif strncmp(model_name,'P948',4)
        dac_mode = 16;
        num_ddr = int32(num_channels/2);
        markers_per_chan = 2;
        max_dac = 65535;
        data_type = 'uint16';
    else
        dac_mode = 16;
        num_ddr = fix(num_channels/2);
        markers_per_chan = 2;
        max_dac = 65535;
        data_type = 'uint16';
    end    
    
    fprintf(1,'DAC waveform format: %d bits-per-point\n',dac_mode); 
    
    half_dac = round(max_dac / 2);
catch ME
    if cType == "LAN"
        inst.Disconnect();
    else
        admin.Close();
    end
    rethrow(ME)
end

awgChan =       2;
dacRes =        dac_mode;
digChan =       2;
numOfFrames =   1;
dualMode =      true; % true = DUAL, false = SINGLE
useDdc =        false;
useDtrTrig =    true;
useExtTrg1 =    false;
trgLevel =      0.2;

awgSRate =      2.5E9;
digSRate =      2.5E9;
lenMult =       2;

digRng =        3; % 1 = LOW, 2 = MED, 3 = HIGH
awgAmp =        0.4;

awgWfm =        1; % 1 = square, 2 = triangle, 3 = cosine, 4 = sine
minCycles =     10;
period =        0.5E-5;
segment =       1;

fprintf(1, 'INITIALIZING SETTINGS\n');

% Reset AWG
inst.SendScpi('*CLS');
inst.SendScpi('*RST');

% Get options using the standard IEEE-488.2 Command
options = pfunc.getOptions(inst);

% Get granularity
granul = pfunc.getGranularity(model_name, options, false);
granul = 64;

% Get Number of Channels
numOfChannels = 2;

% Calculate waveform

fprintf(1, 'Calculating WAVEFORMS\n');

% Sample rate is limited to 2.5GS/s so it is alway possible to generate
% waveforms for all the channels no matter the target Proteus model.
% For the P908X sampling rate can reach 9GS/s in direct mode for all the
% channels. For the P948X models, half the channels are available for
% direct generation over 2.5GS/s without interpolation,
if awgSRate > 2.5E9
    awgSRate = 2.5E9;
end

mySqrWfm = getSquareWfm(    awgSRate,... 
                            minCycles,...
                            period,...
                            granul);

numOfWfms = 4;

awgLength = length(mySqrWfm);

myWfm = zeros(numOfChannels,awgLength);

for wfmNumber = 1:numOfWfms
    myWfm(wfmNumber, :) = mySqrWfm;
   
    % Integration
    mySqrWfm = cumsum(mySqrWfm);
    % DC removal
    mySqrWfm = mySqrWfm - mean(mySqrWfm);
    % Normalization to the -1.0/+1.0 range
    mySqrWfm = mySqrWfm / max(abs(mySqrWfm));
  
end

for wfmNumber = 1:numOfWfms
    myWfm(wfmNumber, (awgLength/2 + 1):awgLength) = 0.0;
end

% SETTING AWG
fprintf(1, 'SETTING AWG\n');

% Set sampling rate for AWG to maximum.
inst.SendScpi([':FREQ:RAST ' num2str(awgSRate)]);
%inst.SendScpi(':TRAC:FORM U8');

%Select Channel
inst.SendScpi(sprintf(':INST:CHAN %d', awgChan));
% DAC Mode set to 'DIRECT" (Default)
inst.SendScpi(':SOUR:MODE DIRECT');

% All segments deleted for current waveform memory bank
inst.SendScpi(':TRAC:DEL:ALL');  

% Waveform Downloading
% *******************

% myWfmSine = 0: (awgLength - 1);
% numOfCycles = awgLength / awgSRate;
% numOfCycles = round(numOfCycles * 6E2);
% 
% mywfmSine = sin(2 * pi * myWfmSine * numOfCycles / awgLength);

fprintf(1, 'DOWNLOADING WAVEFORM FOR CH%d\n', awgChan);

myWfm = myWfm(awgWfm, :);

res = SendWfmToProteus2(    inst,...
                            pfunc,...
                            awgSRate,...
                            awgChan,...
                            segment,...
                            myWfm,...
                            dacRes,...
                            false);

%res = SendWfmToProteus(inst, awgChan, segment, myWfm(awgWfm, :), dacRes);
fprintf(1, 'WAVEFORM DOWNLOADED!\n');
% Select segment for generation
fprintf(1, 'SETTING AWG OUTPUT\n');
inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segment));
% Output volatge set to awgAmp
inst.SendScpi(sprintf(':SOUR:VOLT %d', awgAmp));   
% Activate outpurt and start generation
inst.SendScpi(':OUTP ON');

% The Task Composer is configured to handle a certain number of task
% entries. One task is enough for this application.
inst.SendScpi(sprintf(':TASK:COMP:LENG %d', 1)); 

% Then, each task is defined

% Task to be defined is selected
inst.SendScpi(sprintf(':TASK:COMP:SEL %d', 1));
% The type of task is defined. SINGle is the default so sending this 
% command is not mandatory

inst.SendScpi(':TASK:COMP:TYPE SING');

if useExtTrg1
    inst.SendScpi(':TASK:COMP:ENAB TRG1');
end
    
% The action to take after completing the task is defined. NEXT is the
% default so sending this command is not mandatory
inst.SendScpi(':TASK:COMP:DEST NEXT');

% Assigns segment for task in the sequence        
inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', segment));    
% Next Task is the following task in the table except for the last
% task, which points to task #1    
inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d', 1));    

% Trigger Digitizer  
inst.SendScpi(':TASK:COMP:DTR ON');   

param = 1; 
inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', 1));


% The task table created with the Composer is written to the actual task
% table of teh selected channel
inst.SendScpi(':TASK:COMP:WRIT');
fprintf(1, 'SEQUENCE CREATED!\n');

fprintf(1, 'SETTING AWG OUTPUT\n');
% Select Task Mode for generation    
% Start in task #1 (#1 is the default)    
inst.SendScpi(':FUNC:MODE TASK');
inst.SendScpi(':SOUR:FUNC:MODE:SEGM 1');
% Output volatge set to awgAmp
inst.SendScpi(sprintf(':SOUR:VOLT %d', awgAmp));   
% Activate outpurt and start generation
inst.SendScpi(':OUTP ON');

% TRIGGER SETTING
if useExtTrg1
    inst.SendScpi(':TRIG:SEL TRG1');
    inst.SendScpi(':TRIG:LEV 0.0');
    inst.SendScpi(':TRIG:SLOP POS');
    inst.SendScpi(':TRIG:MODE EVEN');
    inst.SendScpi(':TRIG:STAT ON');
end



%************************************************************
%*                        DIGITIZER                         *
%************************************************************

frameLen = awgLength * digSRate / awgSRate;

frameLen = lenMult * frameLen;
frameLen = ceil(frameLen / 48) * 48;

%frameLen = 270960;

%dtrChan = GetDtrChan(model, awgChan);
dtrChan = awgChan;

if useDtrTrig
else
    dtrChan= 0;
end
holdFlag = false;
while true
    
  % wfmData = GetDigitizerData2( inst,...
%                                 pfunc,...
%                                 cType,...
%                                 dualMode,...
%                                 digChan,...
%                                 digRng,...
%                                 useDdc,...
%                                 0,...
%                                 0.0,...
%                                 1E+9,...
%                                 0.0,...
%                                 1,...
%                                 33888,...
%                                 11363,...
%                                 0); %0); 

    wfmData = GetDigitizerData2( inst,...
                                pfunc,...
                                cType,...
                                dualMode,...
                                digChan,...
                                digRng,...
                                useDdc,...
                                dtrChan,...
                                trgLevel,...
                                digSRate,...
                                0.0,...
                                1,...
                                frameLen,...
                                numOfFrames,...
                                frameLen / 2); %0); %2 * frameLen / 4);
    xData = 0: (frameLen - 1);
    xData = xData / digSRate;
    
    plot(xData, wfmData(1, :));
    if holdFlag
        hold;
        holdFlag = false;
    end
    drawnow;
    
    if numOfFrames > 1
        for i = 2:numOfFrames
            plot(xData, wfmData(i, :));
        end
    end

end

% It is recommended to disconnect from instrument at the end
% It is recommended to disconnect from instrument at the end
if cType == "LAN"
    inst.Disconnect();
else
    admin.CloseInstrument(instId);    
    admin.Close();
end     
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

function result = SendWfmToProteus2(    inst,...
                                        pfunc,...
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
    myWfm = pfunc.myQuantization(myWfm, dacRes, 1);
    
    % Download the binary data to segment   
    prefix = ':TRAC:DATA 0,';

    if (dacRes==16)
        myWfm = uint16(myWfm);
        myWfm = typecast(myWfm, 'uint8');
    else
        myWfm = uint8(myWfm);
    end
    tic;
    %res = inst.WriteBinaryData(':TRAC:DATA ', myWfm);
    res = inst.WriteBinaryData(prefix, myWfm);
    toc
    assert(res.ErrCode == 0);
    
%     if dacRes == 16
%         inst.SendBinaryData(prefix, myWfm, 'uint16');
%     else
%         inst.SendBinaryData(prefix, myWfm, 'uint8');
%     end   
    
    if initialize
        inst.SendScpi(sprintf(':SOUR:FUNC:MODE:SEGM %d', segment));
        % Output voltage set to
        MAX
        inst.SendScpi(':SOUR:VOLT MAX');   
        % Activate outpurt and start generation
        inst.SendScpi(':OUTP ON');        
    end
    
    result = length(myWfm);
end

function wfmData = GetDigitizerData2(   inst,...
                                        pfunc,...
                                        cType,...
                                        dualMode,...
                                        digChan,...
                                        digRng,...
                                        useDdc,...                                   
                                        awgChan,...
                                        trgLevel,...
                                        samplingRateDig,...
                                        cFreq,...
                                        numOfAcq,...
                                        frameLen,...
                                        numberOfFrames,...
                                        preTrig)
    % ----------
    % ADC Config
    % ----------
    
    % Turn on ADC dual-channels mode
    if dualMode
        inst.SendScpi(':DIG:MODE DUAL');
    else
        inst.SendScpi(':DIG:MODE SINGLE');
    end

     % Select digitizer channel:
    if digChan < 1
        digChan = 1;
    elseif digChan > 2
        digChan = 2;
    end

    inst.SendScpi(sprintf(':DIG:CHAN %d', digChan));  

    % Enable acquisition in the digitizer's channels 

    ddcFreq = 0.0;
    if useDdc
        ddcFreq = GetDdcFreq(cFreq, samplingRateDig);
        inst.SendScpi(':DIG:DDC:MODE COMP');   
        %inst.SendScpi(':DIG:DDC:BIND OFF');        
        inst.SendScpi([':DIG:DDC:CFR1 ' num2str(abs(ddcFreq))]);        
    end

    % Set sampling rate
    
    inst.SendScpi(sprintf(':DIG:FREQ %g', samplingRateDig));

    if useDdc        
        %inst.SendCmd(':DIG:DDC:CLKS DIG');
    end       

    if digRng < 1
        digRng = 1;
    elseif digRng > 3
        digRng = 3;
    end

    switch digRng

        case 1
            inst.SendScpi(':DIG:CHAN:RANG LOW'); %MED
        case 2
            inst.SendScpi(':DIG:CHAN:RANG MED'); %MED
        case 3
            inst.SendScpi(':DIG:CHAN:RANG HIGH'); %MED
    end   
        
    
    % Enable acquisition in the digitizer's channels  
    inst.SendScpi(':DIG:CHAN:STATE ENAB');


    actualFrameLen = frameLen;
    actualNumberOfFrames = numberOfFrames;

    if useDdc
        actualFrameLen = 2 * actualFrameLen;
        actualNumberOfFrames = 1 * actualNumberOfFrames;   
    end
    
    % Setup frames layout    
    inst.SendScpi(sprintf(':DIG:ACQ:FRAM:DEF %d, %d',actualNumberOfFrames, actualFrameLen));
    
    % Set channel 1 of the digitizer as it trigger source
    if awgChan == 0   
        inst.SendScpi(sprintf(':DIG:TRIG:SOUR CH%d', digChan)); %digChan
        inst.SendScpi(sprintf(':DIG:TRIG:SELF %f', trgLevel)); %0.025 
    else
        inst.SendScpi(sprintf(':DIG:TRIG:SOURCE TASK%d', awgChan)); %change awgChan
    end
   
    if useDdc
        actualPreTrig = 2 * preTrig;
    else
        actualPreTrig = preTrig;
    end
    actualPreTrig = round( actualPreTrig / 48) * 48;

    inst.SendScpi(sprintf(':DIG:PRET %d', actualPreTrig));

    % Select channel 1
    inst.SendScpi(sprintf(':DIG:CHAN %d', digChan));  

    % Select which frames are filled with captured data 
    %(all frames in this example)
    inst.SendScpi(':DIG:ACQ:FRAM:CAPT:ALL');

    % Delete all wfm memory
    inst.SendScpi(':DIG:ACQ:ZERO:ALL');

    resp = inst.SendScpi(':DIG:DATA:FORM?');
    resp = strtrim(pfunc.netStrToStr(resp.RespStr));

    inst.SendScpi(':DIG:INIT OFF'); 
    
    %tic
    wfmData = zeros(numOfAcq, frameLen);

    for i=1:numOfAcq        
        inst.SendScpi(':DIG:INIT ON');  
        %inst.SendCmd('*TRG');       
        for n = 1:250
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

%         toc

        inst.SendScpi(':DIG:INIT OFF');

        % Define what we want to read 
        % (frames data, frame-header, or both).
        % In this example we read the frames-data
        inst.SendScpi(':DIG:DATA:TYPE FRAM');
        inst.SendScpi(':DIG:DATA:SEL ALL');

        % Read binary block
        resp = inst.SendScpi(':DIG:DATA:SIZE?');
        resp = strtrim(pfunc.netStrToStr(resp.RespStr));
        num_bytes = str2double(resp);        
        
        if useDdc
            adc_data_type = 'uint32';
            fprintf(1, 'WFM Transfer from Proteus Starts\n');
            
            samples = inst.ReadBinaryData(':DIG:DATA:READ?', 'uint32');
            fprintf(1, 'WFM Transfer from Proteus Finished. Time:\n');
            fprintf(1, 'WFM Transfer Time:\n');
%             toc;   
            samples = int32(samples) - 16384;
            samples = double(samples(1:2:length(samples)))  + 1i * double(samples(2:2:length(samples)));
        else
            adc_data_type = 'uint16';
            if cType == "LAN"
                %tic;
                samples = inst.ReadBinaryData(':DIG:DATA:READ?', adc_data_type);
                samples = cast(samples,adc_data_type);
                %toc
            else
                % because read format is UINT16 we divide byte number by 2
                wavlen = floor(num_bytes / 2);
            
                % allocate NET array
                netArray = NET.createArray('System.UInt16', wavlen);%wavelen
                % read the captured frame
                %tic
                res = inst.ReadMultipleAdcFrames(digChan - 1, 1, numberOfFrames, netArray); % 0,1, numberofframes, netarray
                %toc
                assert(res == 0);
                
                % cast to matlab vector
    	        samples = uint16(netArray);
                                
                % deallocate the NET array
                delete(netArray);
            end            
            samples = int16(samples) - 2048;
        end

        for j=1:numberOfFrames
            wfmData((i - 1) * numberOfFrames + j,:) = samples(((j-1) * frameLen + 1):(j * frameLen));
        end        
    end
end

function ddcFreq = GetDdcFreq(carrierF, adcSr)
    ddcFreq = carrierF;
    while ddcFreq > adcSr/2
        ddcFreq = ddcFreq - adcSr;
    end
end

function dtrChan = GetDtrChan(model, outChNum)

    if outChNum > 4
        outChNum = 4;
    elseif outChNum < 1
        outChNum = 1;
    end
    
    dtrChan = outChNum;

    model = upper(model);
    model = strtrim(model);
    model = split(model, ',');

    if contains(model(1), "TABOR")
        model = model(2);
    else
        model = model(1);
    end

    if contains(model, 'D')
        switch outChNum
            case 1
                dtrChan = 3;

            case 2
                dtrChan = 4;

            case 3
                dtrChan = 1;

            case 4
                dtrChan = 2;
        end
    end
end