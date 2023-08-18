% This version runs a pulse-acquire FID
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
    
%     Loop to repeatedly wait for messages and send replies
%     Break or Ctrl+C to get out of loop
    while ( connect==on )
        % Wait for message
        while(u2.BytesAvailable == 0)
            % If no bytes in u2 buffer, wait 10ms then check again
            pause(0.01);
        end
        %         cmdBytes = fread(u2);
        readBytes = fscanf(u2);
        dataBytes=1;counter=1;
        while ~isempty(readBytes)
            [tempBytes,readBytes]=strtok(readBytes,',');
            cmdBytes(counter)=str2num(tempBytes);
            counter=counter+1;
        end
        %dataBytes=strtok(dataBytes,',');
        cmdByte=cmdBytes(1);
        
        % 1 - Init
        % 2 - Aquire on Trig
        % 3 - Measure
        % 4 - Cleanup - prepare for next cycle
        % 5 - Disconnect Intrument
        % 6 - Program MW chirp waveform
        % 7 - Play MW chirp waveform
        
        %measure_time=str2double(dataBytes); %in sec
        switch cmdByte
            case 1 % Init
                
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
    
    %inst.SendScpi(sprintf(':DIG:CHAN %d', 1)); 
    %inst.SendScpi(':DIG:ACQ:FREE');    
    %inst.SendScpi(':DIG:DDC:CLKS AWG');
    %inst.SendScpi(':DIG:DDC:MODE COMP');T
    %inst.SendScpi(':DIG:DDC:CFR1 0');
    %inst.SendScpi(':DIG:DDC:PHAS1 90');
    %inst.SendScpi(':DIG:CHAN:RANG HIGH');
    % Enable acquisition in the digitizer's channels  
    %inst.SendScpi(':DIG:CHAN:STATE ENAB'); 
    
    fprintf('ADC Configured\n');
    fprintf('Clocks synced\n');
                
            case 2 % Aquire on trig
                
    % ---------------------------------------------------------------------
    % RF Pulse Config
    % ---------------------------------------------------------------------
    
%     pulse_name = ['init_pul'];
    amps = [0.5];
    frequencies = [0];
    lengths = [100e-6];
    phases = [0];
    mods = [0]; %0 = square, 1=gauss, 2=sech, 3=hermite 
    spacings = [50000e-6];
    markers = [1]; %always keep these on
    markers2 = [0];
    trigs = [1]; %acquire on every "pi" pulse
    
    reps = [1];
    repeatSeq = [1]; % how many times to repeat the block of pulses
    
                pw = cmdBytes(2)*1e-6;
                lengths(1) = pw;
%                 tof = -1000*cmdBytes(2);
                tof = -1000*(25.0613);
                
                ch=1;
                initializeAWG(ch);
                clearPulseDict();
                clearBlockDict();
                
                defPulse('init_pul', amps(1), mods(1), lengths(1), phases(1), spacings(1));
                defBlock('FID', {'init_pul'}, reps(1), markers(1), trigs(1));
                makeBlocks({'FID'}, ch, repeatSeq);
                    
                setNCO_IQ(ch, 75.38e6+tof, 0);
                inst.SendScpi(sprintf(':DIG:DDC:CFR2 %d', 75.38e6+tof));
                
% %                 need to modify Marker #2 to show up during acquisition
%                 % ## final segment which will contain the marker
%                 chNum = 1;
%                 segNum = 4;
%                 
%                 % ## how long the "on" portion needs to be
%                 onLength = 10e-6;
%                 
%                 % ## when does the marker start after the last pulse
%                 buffer = 10e-6;
%                 
%                 % ## select the final segment
%                 cmd = sprintf(':INST:CHAN %d',chNum);
%                 inst.SendScpi(cmd);
%                 cmd = sprintf(':TRAC:SEL %d',segNum);
%                 inst.SendScpi(cmd);
%                 
%                 % ## get the length of the segment
%                 query = inst.SendScpi(':TRAC:DEF?');
%                 segLen = pfunc.netStrToStr(query.RespStr);
%                 mkrsegment_length = str2num(segLen);
% %                 mkrsegment_length = floor(mkrsegment_length/4);
%                 
%                 % ## make a new segment
%                 mkr_vector_2 = zeros(mkrsegment_length,1);
%                 mkr_vector_1 = zeros(mkrsegment_length,1);
%                 onLength_points = floor(onLength*sampleRateDAC/32)*8;
%                 buffer_start = floor(buffer*sampleRateDAC/32)*8;
%                 
%                 % ## make the marker
% %                 mkr_vector_on = ones(onLength_points,1);
%                 mkr_vector_2((buffer_start+1) : buffer_start+onLength_points) = 1;% mkr_vector_on;
%                 % proteus.inst.timeout = 30000
%                 
%                 % # Send the binary-data with *OPC? added to the beginning of its prefix.
%                 mkr_vector = mkr_vector_1 + 2*mkr_vector_2;
% %                 mkr_vector = mkr_vector(1:2:length(mkr_vector)) + 16 * mkr_vector(2:2:length(mkr_vector));
%                 mkr_vector = uint8(mkr_vector);
%                 % inst.WriteBinaryData('*OPC?; :MARK:DATA', mkr_vector);
%                 inst.WriteBinaryData(':MARK:DATA 0,', mkr_vector);
%                 
% %                 % % # Set normal timeout
% %                 % proteus.inst.timeout = 10000
% %                 cmd = ':MARK:SEL 1';
% %                 inst.SendScpi(cmd);
% %                 cmd = ':MARK:STAT ON';
% %                 inst.SendScpi(cmd);
% %                 % proteus.checkForError()

                
                fprintf('Calculate and set data structures...\n');
                
                numberOfPulses_total = 1;

                Tmax=1; % will be 1 for FID
                
                tacq=50000;

                numberOfPulses=1; %in 1 second %will be 1 for FID
                loops=Tmax;
                    
                readLen = round2((tacq+2)*1e-6*2.7/(16*1e-9),96)-96;
                
                offLen = 0;
                rc = inst.SendScpi(sprintf(':DIG:ACQ:DEF %d, %d',numberOfPulses*loops, 2 * readLen));
                assert(rc.ErrCode == 0);
                
                inst.SendScpi(sprintf(':DIG:CHAN %d', adcChanInd))
                %rc = inst.SendScpi(':DIG:TRIG:SOUR TASK1'); %digChan
                %assert(rc.ErrCode == 0);
                %rc = inst.SendScpi(sprintf(':DIG:TRIG:SELF %f', 0.025)); %0.025 
                %assert(rc.ErrCode == 0);
                %rc = inst.SendScpi(sprintf(':DIG:TRIG:AWG:TDEL %f', lengths(2) + 8e-6));
                %assert(rc.ErrCode == 0);
                
                %inst.SendScpi(sprintf(':DIG:CHAN 1'))
                %rc = inst.SendScpi(':DIG:TRIG:SOUR TASK1'); %digChan
                %assert(rc.ErrCode == 0);
%                 rc = inst.SendScpi(sprintf(':DIG:TRIG:AWG:TDEL %f', 0));
%                 assert(rc.ErrCode == 0);
                
                rc = inst.SendScpi(':DIG:TRIG:SOUR EXT'); %digChan
                assert(rc.ErrCode == 0);
                %inst.SendScpi(':DIG:TRIG:TYPE GATE');
                rc = inst.SendScpi(':DIG:TRIG:SLOP NEG');
                assert(rc.ErrCode == 0)
                rc = inst.SendScpi(':DIG:TRIG:LEV1 0.01');
                assert(rc.ErrCode == 0)
                rc = inst.SendScpi(sprintf(':DIG:TRIG:DEL:EXT %f', 6e-6)); % external trigger delay
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
                
                
            case 3 % Measure
                
                %inst.SendScpi(':DIG:INIT ON');
                inst.SendScpi(sprintf(':DIG:CHAN 2'))
                fprintf('Triggering pulse sequence\n');
                rc = inst.SendScpi('*TRG');
                assert(rc.ErrCode == 0);
                
                n=0;
                
               % pause(Tmax+3);
                
                for n = 1:700
                    
                    resp = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
                    resp = strtrim(pfunc.netStrToStr(resp.RespStr));
                     items = split(resp, ',');
                     items = str2double(items);
                     if length(items) >= 3 && items(2) == 1
                         break
                     end
                     if mod(n, 10) == 0
                        fprintf('%d. %s Time:\n', fix(n / 10), resp);
%         %                 toc                
                     end
                     pause(0.1);
                    
                 end

                rc = inst.SendScpi(':DIG:INIT OFF');
                assert(rc.ErrCode == 0);
%                 rc = inst.SetAdcCaptureEnable(on);
%                 assert(rc == 0);
                
%                 rc = inst.ReadAdcCaptureStatus();
                
                %                 % Wait until the capture completes
                %                 status = inst.ReadAdcCompleteFramesCount();
                %                 while status ~= numberOfPulses_total
                %                     pause(0.01);
                %                     status = inst.ReadAdcCompleteFramesCount();
                %                 end
                
                %                 readSize = uint64(readLen);
                %                 readOffset = uint64(offLen);
                
                chanIndex = adcChanInd - 1;
                pulseAmp = [];
                relPhase = [];
                
                netArray = NET.createArray('System.UInt16', 4*readLen*numberOfPulses); %total array -- all memory needed
                %netArray2 = NET.createArray('System.UInt16', 4*readLen*numberOfPulses);
                
                power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
                padded_len= 2^(power_of_2) ;%2^15;
                %padded_len = readLen;
                %padded_len = 4096;
                dF = sampleRate/16/padded_len; %set the discretization of freq in terms of sampleRate
                f = -sampleRate/16/2:dF:sampleRate/16/2-dF; %sampleRate sets the 'bandwidth'
                
                %[v,b1]=min(abs(f-20.00706e6)); %picks out the 20MHz component (Varian)
                %[v,b2]=min(abs(f+20.00706e6));
                
                [v,b1]=min(abs(f-75.35e6)); %picks out the 75 MHz component (Tecmag)
                [v,b2]=min(abs(f+75.35e6));
                
                eff_read=100:readLen-100;
                cyclesPoints = 50;
                fprintf('Shuttle complete\n')
                fprintf('Transfering aquired data to computer....\n')
                for n = 1:loops
                    fprintf('Start Read %d .... ', n);
                    firstIndex = ((n-1)*numberOfPulses)+1;
                    
                    % Define what we want to read 
                    % (frames data, frame-header, or both).
                    % In this example we read the frames-data
                     rc = inst.SendScpi(':DIG:DATA:TYPE FRAM');
                     assert(rc.ErrCode == 0);
                     rc = inst.SendScpi(':DIG:DATA:SEL ALL');
                     assert(rc.ErrCode == 0);

                    % Read binary block
                    resp = inst.SendScpi(':DIG:DATA:SIZE?');
                    resp = strtrim(pfunc.netStrToStr(resp.RespStr));
                    num_bytes = str2double(resp);
                    
                    wavlen = floor(num_bytes / 2);
%                     netArray = NET.createArray('System.UInt16', wavlen);
                    nFrames = 100;
%                     netArray2 = NET.createArray('System.UInt16', nFrames * readLen);
                    
                    tic             
%                     
                                                       
                    rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, numberOfPulses, netArray); %this is where the device reads
                    %rc = inst.ReadMultipleAdcFrames(0, firstIndex, numberOfPulses, netArray2); %this is where the device reads
                    
                    assert(rc == 0);
                    samples = uint16(netArray); %get the data (1s chunk)
                    samples = samples(1:2:length(samples));
                    samples = int16(samples)-16384;
                    wfmLength = length(samples);
                    samplesI = double(samples(1:2:wfmLength));
                    samplesQ = double(samples(2:2:wfmLength));
                    samples = samplesI + 1.0i * samplesQ;
                    %samples2 = samples2(1:2:length(samples2));
                    %samples2 = samples2(1:2:length(samples2));
                    fprintf('Read %d Complete\n', n);
                    toc
                    
                    tic
                    %delete mem
                    fprintf('Clear mem %d .... ', n);
                    
%                     rc =inst.WipeMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, 0);
%                     assert(rc == 0);
                    fprintf('Clear mem %d Complete\n', n);
                    toc
                    
                    tic
                    fprintf('Starting iteration %d data processing....', n);
                    pulses = reshape(samples, [], numberOfPulses); % reshape samples into a more useful array (2 dimensions)
                    %pulses2 = reshape(samples2, [], numberOfPulses);
                    
                    clear samples;
                    clear samples2;
%                     tWin = flattopwin()
                    for i = 1:numberOfPulses
                        pulse = pulses(:, i);
                        
                        if n == 1
                            if i == 1
                                figure(6);clf;
                                plot(pulse);
                                figure(7);clf;
                                plot(f,abs(fftshift(fft(pulse,padded_len))));
                                hold on;
                                yline(2048);
                            end
                        end

                        idx = i+(numberOfPulses*(n-1));
                        realMean = mean(real(pulse));
                        imagMean = mean(imag(pulse));
                        pulseAmp(idx) = abs(realMean + 1.0i*imagMean);
                        relPhase(idx) = angle(realMean + 1.0i*imagMean);
                    end

                    clear pulses;
                    fprintf('Data processing iteration %d complete!\n', n);
                    toc
                end
                
                time_vec = 1:length(pulse);
                time_axis = time_vec./sampleRate;
                complex_acq = pulse;
%                 
                %fn=dataBytes; %filename
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                % Save data
                fprintf('Writing data to Z:.....\n');
                save(['Z:\' fn],'complex_acq','time_axis');
                fprintf('Save complete\n');
                
            case 4 % Cleanup, save and prepare for next experiment
                rc = inst.SendScpi(':DIG:INIT OFF');
                assert(rc.ErrCode == 0);
%                 rc = inst.SetAdcCaptureEnable(off);

                % Free the memory space that was allocated for ADC capture
                % Delete all wfm memory
                rc = inst.SendScpi(':DIG:ACQ:FREE');
                assert(rc.ErrCode == 0);
%                 rc = inst.FreeAdcReservedSpace();
%                 assert(rc == 0);
                
                %                 clear pulseAmp time_axis;
                fprintf('Complete ready for next aquisition!\n');
                
            case 5 % Shutdown disconnect instr
                
                connect = 0;
            
            case 6 % Program MW chirp waveform
                
                % ---------------------------------------------------------------------
                % Download chirp waveform of 1024 points to segment 1 of channel 1
                % ---------------------------------------------------------------------                
                global sampleRateDAC
                sampleRateDAC = 9e9
                awg_center_freq = cmdBytes(2);
                awg_bw_freq = cmdBytes(3);
                awg_amp = cmdBytes(4);
                sweep_freq = cmdBytes(5);
                sweep_sigma = cmdBytes(6);
                symm = cmdBytes(7);
%                 srs_freq = cmdBytes(8); % should be 0.3625e9
                srs_freq = 0.3625e9; % new value for good chirp
                srs_amp = cmdBytes(9);
                pol_times = [cmdBytes(10) cmdBytes(11) cmdBytes(12) cmdBytes(13) cmdBytes(14) cmdBytes(15)];
                pol_times = nonzeros(pol_times);
                starting_pol_sign = cmdBytes(16);
                
                if starting_pol_sign == 1
                    starting_pol = '+';
                elseif starting_pol_sign ==-1
                    starting_pol = '-';
                else
                    disp('starting polarization value is wrong (only +1 or -1)')
                    break
                end
                
       
%                 awg_center_freq = 3.775e9;
%                 awg_bw_freq = 24e6;
%                 awg_amp = 1.2;
%                 sweep_freq = 750;
%                 sweep_sigma = 0.1;
%                 symm = 0;
%                 srs_freq = 0.3625e9; % new value for good chirp
%                 srs_amp = 3;
%                 pol_times = [t1 t2 t3];
%                 pol_times = nonzeros(pol_times);
%                 starting_pol_sign = 1;
%                 
                inst.SendScpi("*CLS")
                inst.SendScpi("*RST")
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                % check if its the good variable bc im not sure
                inst.SendScpi([':FREQ:RAST ' num2str(sampleRateDAC)]);
                assert(res.ErrCode == 0);
                
                inst.SendScpi(":INIT:CONT OFF");
                assert(res.ErrCode == 0);
                
                inst.SendScpi(":INIT:CONT ON");
                assert(res.ErrCode == 0);
                                
                rampTime = 1/sweep_freq;
                fCenter = awg_center_freq - srs_freq;
                fStart = fCenter - 0.5*awg_bw_freq;
                disp(['fstart = ' num2str(fStart)]);
                fStop = fCenter + 0.5*awg_bw_freq;
                dt = 1/sampleRateDAC;
                
                chirps{1}.dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
                chirps{2}.dacSignal = fliplr(chirps{1}.dacSignal);
                chirps{1}.segLen = length(chirps{1}.dacSignal);
                chirps{2}.segLen = length(chirps{2}.dacSignal);
                chirps{1}.segm = 1;
                chirps{2}.segm = 2;
                chirps{1}.srs_freq = 0.3625e9;
                chirps{2}.srs_freq = 0.3625e9;
                fprintf('waveform length - ');
                fprintf(num2str(length(chirps{1}.dacSignal)));
                fprintf('\n') ;


task_list = build_tasktable(inst,pol_times,chirps,sampleRateDAC,'first sign',starting_pol);

% % Play seg 1
% res = inst.SendScpi(':SOUR:FUNC:MODE:SEGM 1');
% assert(res.ErrCode == 0);
 
res = inst.SendScpi(':SOUR:VOLT MAX');
assert(res.ErrCode == 0);

res = inst.SendScpi(':OUTP ON');
assert(res.ErrCode == 0);

inst.SendScpi(':TRIG:ACTIVE:SEL TRG2');
inst.SendScpi(':TRIG:LEV 1.0');
inst.SendScpi(':TRIG:ACTIVE:STAT ON');

% create the task table
create_task_table(inst,task_list, sampleRateDAC);

% write task table
inst.SendScpi(':TASK:COMP:WRITE');
fprintf(1, 'SEQUENCE CREATED!\n');

fprintf(1, 'SETTING AWG OUTPUT\n');

inst.SendScpi(':SOUR:FUNC:MODE TASK');

res = inst.SendScpi(':OUTP ON');
assert(res.ErrCode == 0);
                
            case 7 % Play MW chirp waveform
                
                % ---------------------------------------------------------------------
                % Play segment 1 in channel 1
                % ---------------------------------------------------------------------
                
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':SOUR:FUNC:MODE TASK');
                assert(res.ErrCode == 0);

%               inst.SendScpi(':INIT:CONT ON');
%               assert(res.ErrCode == 0);
%                 res = inst.SendScpi(':SOUR:FUNC:MODE ARB');
%                 assert(res.ErrCode == 0);
%                 res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 1');
%                 assert(res.ErrCode == 0);
                
                amp_str = ':SOUR:VOLT MAX';
                res = inst.SendScpi(amp_str);
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':OUTP ON');
                assert(res.ErrCode == 0);
                
                fprintf('Waveform generated and playing\n');
                
            case 8
                
                % Disable MW chirp output
                res = inst.SendScpi(':OUTP OFF');
                assert(res.ErrCode == 0);
                
                fprintf('MW Chirp Waveform stopped playing (on purpose)\n');
                
            case 9
                
                % ---------------------------------------------------------------------
                % Play segment 2 in channel 1
                % ---------------------------------------------------------------------
                
%                 res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
%                 assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 2');
                assert(res.ErrCode == 0);
                
                %amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
                amp_str = [':SOUR:VOLT MAX'];
                res = inst.SendScpi(amp_str);
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':OUTP ON');
                assert(res.ErrCode == 0);
                
                fprintf('Waveform in Segment 2 playing\n');
                
        end % end switch
        
    end % end while
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd - server stopped!! \nInstrument Error Status = %s\n', netStrToStr(res.RespStr));
    
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
 