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
    fprintf('initializing Tektronix AFG 31000\n');
    tek = Tektronix_AFG_31000("USB0::0x0699::0x0355::C019986::INSTR");
    fprintf("Tektronix Initialization complete\n");
    
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
    tek.output_off();
                
            case 2 % Aquire on trig
                
    % ---------------------------------------------------------------------
    % RF Pulse Config
    % ---------------------------------------------------------------------
    sampleRateDAC_freq = 675000000;
    granularity = 64;
    
    %% DEFINE PULSE LENGTH
    pi = cmdBytes(3)*1e-6;
    index = cmdBytes(2);
    
    %% DEFINE PULSE SEQUENCE PARAMETERS
    amps = [1 1 1 1];
    frequencies = [0 0 0 0];
    lengths = [pi/2 pi pi/2 pi/2];
    num_init_spin_lock = 6000;
    lengths = round_to_DAC_freq(lengths, sampleRateDAC_freq, 64);
    phases = [0 0 90 90];
    mods = [0 0 0 0]; %0 = square, 1=gauss, 2=sech, 3=hermite
    
    spacings = [5e-6 36e-6 36e-6 5e-6];
    spacings = round_to_DAC_freq(spacings, sampleRateDAC_freq, 64);
    markers = [1 1 1 1]; %always keep these on => turns on the amplifier for the pulse sequence
    trigs = [0 1 1 0];
    
    T_arr = round_to_DAC_freq([2.0e-3, 2.4e-3], sampleRateDAC_freq, granularity);
    [T1, T2] = deal(T_arr(1), T_arr(2));
    
    seq_time = 10;
    % rounding happens in the "get_DTQCs_seq" function
    seg_idx_l = get_DTQCs_seq(sampleRateDAC_freq, granularity, lengths, spacings, T1, T2, seq_time);
    
    fprintf("setting up pulse blaster sequence\n");
    PB = containers.Map('KeyType', 'double', 'ValueType', 'any');
    ch3 = 3;
    PB_seg1 = zeros(2, 2);
    start_time = lengths(1) + spacings(1) + ...
               (lengths(3) + spacings(3))*num_init_spin_lock + lengths(2) + spacings(2) + lengths(2)/2 + 100*T2;
    
    % Below means that PB outputs 0V for duration start_time
    % and PB outputs 2.7V (TTL +) for 100e-6
    [PB_seg1(1,1), PB_seg1(2,1)] = deal(0, 1);
    [PB_seg1(1,2), PB_seg1(2,2)] = deal(start_time, 100e-6);
    fprintf(sprintf("This is AC start time: %d \n", start_time));
    PB(ch3) = PB_seg1;
    initializeAWG(ch3);
    fprintf("downloading pulseblaster sequence \n");
    generate_PB(PB, sampleRateDAC, inst);
    fprintf("PB download finished \n");
    setNCO_IQ(ch3, 0, 0);
    
    freq = 1/T2;
    AC_dict.freq = freq;
    AC_dict.DC_offset = 0;
    AC_dict.Vpp = 0.5;
    AC_dict.phase = -90;
    fprintf(sprintf("This is AC frequency: %d \n", AC_dict.freq));
    fprintf(sprintf("This AC Vpp voltage: %d \n", AC_dict.Vpp));
    fprintf(sprintf("This is AC DC offset: %d \n", AC_dict.DC_offset));
    fprintf(sprintf("This AC phase: %d \n", AC_dict.phase));
    
    
    % The number of readout windows 2nd and 3rd segment are only ones that
    % have readout window.
    numberOfPulses_total = num_init_spin_lock;
    for i = (1: length(seg_idx_l))
        if seg_idx_l(i, 1) == 2 || seg_idx_l(i, 1) == 3
           numberOfPulses_total = numberOfPulses_total + seg_idx_l(i, 2);
        end
    end
    fprintf("This is the number of total readout window: %d \n", numberOfPulses_total);
    
                tof = cmdBytes(6);
                
                ch=1;
                initializeAWG(ch);
                clearPulseDict();
                clearBlockDict();
                
                generateDTQCseqIQ(ch, amps, frequencies, lengths, phases, spacings, markers, trigs, seg_idx_l, num_init_spin_lock);
                assert(sampleRateDAC_freq == sampleRateDAC, "The two sampleRateDAC frequency should be the same");
                    
                setNCO_IQ(ch, 75.38e6+tof, 0);
                fprintf("snyching Tabor's PB and Pseq \n");
                
                inst.SendScpi(sprintf(':INST:CHAN %d',ch));
                inst.SendScpi(':TRIG:COUPLE ON');
                inst.SendScpi(':TRIG:CPU:MODE LOCAL');
                inst.SendScpi(':TRIG:SOUR:ENAB CPU');
                inst.SendScpi(':TRIG:SEL CPU');
                inst.SendScpi(':TRIG:STAT ON');
                resp = inst.SendScpi(':SYST:ERR?');
                
                inst.SendScpi(sprintf(':DIG:DDC:CFR2 %d', 75.38e6+tof));
                
                fprintf('Calculate and set data structures...\n');

                
                Tmax=cmdBytes(4);
                
                %fix the window as 2us.
                tacq=cmdBytes(5);
                
                fprintf("This is tacq: %d \n", tacq);
%                 tacq=128;
%                 tacq=64;
%                 tacq=96;
                
                numberOfPulses= floor(numberOfPulses_total/Tmax); %in 1 second %will be 1 for FID
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
                rc = inst.SendScpi(':DIG:TRIG:LEV1 1.0');
                assert(rc.ErrCode == 0)
                rc = inst.SendScpi(sprintf(':DIG:TRIG:DEL:EXT %f', 12e-6)); % external trigger delay
                assert(rc.ErrCode == 0)
                
                fprintf('Instr setup complete and ready to aquire\n');
                
                %netArray = NET.createArray('System.UInt16', 2* readLen*numberOfPulses); %total array -- all memory needed
                
%                 rc = inst.SetAdcAcquisitionEn(on,off);
%                 assert(rc == 0);

                % Setup frames layout    

                
%                  rc = inst.SendScpi(':DIG:ACQ:CAPT 1,-1');   
%                  assert(rc.ErrCode == 0);
%                  rc = inst.SendScpi(':DIG:ACQ:ZERO 1,-1');
%                  assert(rc.ErrCode == 0);
                 rc = inst.SendScpi(':DIG:ACQ:FRAM:CAPT:ALL');   
                 assert(rc.ErrCode == 0);
                 rc = inst.SendScpi(':DIG:ACQ:ZERO:ALL');
                 assert(rc.ErrCode == 0);

%                 rc = inst.SetAdcFramesLayotrigut(numberOfPulses*loops, readLen); %set memory of the AWG
    %                 assert(rc == 0);
%                  rc = resp = inst.SendScpi(':DIG:DATA:FORM?');
%                  assert(rc == 0);
%                  resp = strtrim(pfunc.netStrToStr(resp.RespStr));
    
                fprintf('Waiting... Listen for Shuttle\n');
                rc = inst.SendScpi(':DIG:INIT OFF'); 
                assert(rc.ErrCode == 0);
                rc = inst.SendScpi(':DIG:INIT ON');
                assert(rc.ErrCode == 0);
                
                fprintf("set Tektronix 31000 as burst mode \n");
                ncycles = min(1e6, round(60*AC_dict.freq));
                if AC_dict.Vpp~=0 || AC_dict.DC_offset~=0
                    tek.burst_mode_trig_sinwave(AC_dict.freq, AC_dict.Vpp,...
                        AC_dict.DC_offset, AC_dict.phase, ncycles);
                end
                fprintf("setting done\n");
               
%                 rc = inst.SetAdcCaptureEnable(on);
%                 assert(rc == 0);

%                 resp1 = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
%                 resp1 = strtrim(pfunc.netStrToStr(resp1.RespStr));
%                 pause(0.1);
%                 resp2 = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
%                 resp2 = strtrim(pfunc.netStrToStr(resp2.RespStr));
%                 pause(0.1);
%                 resp3 = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
%                 resp3 = strtrim(pfunc.netStrToStr(resp3.RespStr));
                
                
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
                    
                    if savealldata
                        pulsechunk = int16(pulses);
                        a = datestr(now,'yyyy-mm-dd-HHMMSS');
                        fn = sprintf([a,'_Proteus','_chunk', num2str(n)]);
                        % Save data
                        fprintf('Writing data to Z:.....\n');
                        save(['Z:\' fn],'pulsechunk');
                        %writematrix(pulses);
                    end
                    
                    if savesinglechunk
                        if  n==1 %determines which chunk will be saved
                            pulsechunk = int16(pulses);
                            a = datestr(now,'yyyy-mm-dd-HHMMSS');
                            fn = sprintf([a,'_Proteus_chunk', num2str(n)]);
                            % Save data
                            fprintf('Writing data to Z:.....\n');
                            save(['Z:\' fn],'pulsechunk');
                            %writematrix(pulses);
                        end
                    end
                    
                    clear samples;
                    clear samples2;
%                     tWin = flattopwin()
                    for i = 1:numberOfPulses
                        pulse = pulses(:, i);
                       
                        
                        %                         pulse = pulse(1024:readLen);
                        %                         readLen=length(pulse);
                        if savemultwind %save multiple consecutive windows of raw data if true
                            if n==2 & i<=8
                                pulsechunk = int16(pulse);
                                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                                fn = sprintf([a,'_Proteus_chunk', num2str(n) '_' num2str(i)]);
                                % Save data
                                fprintf('Writing data to Z:.....\n');
                                save(['Z:\' fn],'pulsechunk');
                            end
                        end
                        
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
                        if n == 1
                            if i == 2
                                figure(8);clf;
                                plot(pulse);
                                figure(9);clf;
                                plot(f,abs(fftshift(fft(pulse-mean(pulse),padded_len))));
                                hold on;
                                yline(2048);
                            end
                        end
                        if n == 58
                            if i == 9708
                                figure(10);clf;
                                plot(pulse);
                                figure(11);clf;
                                plot(f,abs(fftshift(fft(pulse-mean(pulse),padded_len))));
                                hold on;
                                yline(2048);
                            end
                        end

                        idx = i+(numberOfPulses*(n-1));
                        realMean = mean(real(pulse));
                        imagMean = mean(imag(pulse));
                        pulseAmp(idx) = abs(realMean + 1.0i*imagMean);
                        relPhase(idx) = angle(realMean + 1.0i*imagMean);
                        
                        
%                         pulseDC = pulse - pulseMean; % remove DC
%                         
%                         X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
%                         linFFT = (abs(X)/readLen);
%                         %                         [amp, loc] = max(linFFT);
%                         amp=(linFFT(b1) + linFFT(b2));
%                         phase1=angle(X(b1));
%                         phase2=angle(X(b2));
%                         phase=phase1-phase2;
%                         %             linFFT_vec(:,i)=linFFT;
%                         idx = i+(numberOfPulses*(n-1));
%                         pulseAmp(idx) = amp;
% %                         relPhase(idx) = phase;
%                         relPhase(idx) = phase2;
                    end
%                     if n == 1
%                         figure(8);clf;
%                         plot(pulse);
%                         hold on;
%                         yline(2048);
%                     end
                    clear pulses;
                    fprintf('Data processing iteration %d complete!\n', n);
                    toc
                end
                
                %creating time_axis --- need to account for all the pulses
                t_unit = granularity * vpa(1/sampleRateDAC_freq);
                fprintf("generating time axis... \n");
                time_axis = zeros(1, numberOfPulses_total);
                time_idx = 1;
                curr_t = 0;
                for idx1 = (1:num_init_spin_lock)
                    time_axis(time_idx) = idx1*(lengths(3) + spacings(3));
                    time_idx = time_idx + 1;
                end
                curr_t = num_init_spin_lock * (lengths(3) + spacings(3));
                for idx2 = (1: length(seg_idx_l))
                    if seg_idx_l(idx2, 1) == 2
                        curr_t = curr_t + lengths(2) + spacings(2);
                        time_axis(time_idx) = curr_t;
                        time_idx = time_idx + 1;
                    elseif seg_idx_l(idx2, 1) == 3
                        num_x_pulse = seg_idx_l(idx2, 2);
                        time_axis(time_idx: time_idx + num_x_pulse - 1) = curr_t + (lengths(3) + spacings(3))*(1:num_x_pulse);
                        time_idx = time_idx + num_x_pulse;
                        curr_t = curr_t + (lengths(3) + spacings(3)) * num_x_pulse;
                    elseif seg_idx_l(idx2, 1) == 4
                        curr_t = curr_t + (lengths(4) + spacings(4))*seg_idx_l(idx2, 2);
                    elseif seg_idx_l(idx2, 1) == 5
                        curr_t = curr_t + t_unit*seg_idx_l(idx2, 2);
                    end
                end
                fprintf("Finished generating time axis... \n");
                  %drop first point -- NOT ANYMORE
%                 time_axis(1)=[];pulseAmp(1)=[];relPhase(1)=[];
                phase_base = mean(relPhase(1000:2000)); % take average phase during initial spin-locking to be x-axis
                relPhase = relPhase - phase_base; % shift these values so phase starts at 0 (x-axis)
                try
                    start_fig(12,[5 1]);
                    p1=plot_preliminaries(time_axis,(relPhase),2,'noline');
                    set(p1,'markersize',1);
                    plot_labels('Time [s]', 'Phase [au]');
                    
%                     start_fig(1,[3 2]);
%                     p1=plot_preliminaries(time_axis,pulseAmp,1,'nomarker');
%                     set(p1,'linewidth',1);
%                     set(gca,'ylim',[0,max(pulseAmp)*1.05]);
%                     set(gca,'xlim',[0,25e-3]);
%                     plot_labels('Time [s]', 'Signal [au]');
                    
                    start_fig(1,[5 2]);
                    p1=plot_preliminaries(time_axis,pulseAmp,1,'noline');
                    set(p1,'markersize',1);
                    set(gca,'ylim',[0,max(pulseAmp)*1.05]);
                    plot_labels('Time [s]', 'Signal [au]');
                    
                    start_fig(2,[5 2]);
                    p1=plot_preliminaries(time_axis,zeros(1,length(time_axis)),5,'nomarker');
                    set(p1,'linestyle','--'); set(p1,'linewidth',1);
                    p1=plot_preliminaries(time_axis,pulseAmp.*cos(relPhase),1,'noline');
                    set(p1,'markersize',1);
                    set(gca,'ylim',[-max(pulseAmp)*1.05,max(pulseAmp)*1.05]);
                    plot_labels('Time [s]', 'Signal [au]');
                    
%                     start_fig(13,[5 1]);
%                     p1=plot_preliminaries(time_axis,pulseAmp,1,'nomarker');
%                     set(p1,'linewidth',0.5);
%                     set(gca,'xlim',[5-8e-3,5+30e-3]);
%                     plot_labels('Time [s]', 'Signal [au]');

                catch
                    disp('Plot error occured');
                end
                
                %fn=dataBytes; %filename
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                % Save data
                fprintf('Writing data to Z:.....\n');
                save(['Z:\' fn],'pulseAmp','time_axis','relPhase','AC_dict','lengths',...
                    'phases','spacings','num_init_spin_lock','trigs','start_time');
                fprintf('Save complete\n');
                tek.output_off()
                
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
for iter = (1:10)
    Pines_write(2021, '6');
end 
                
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
                for iter = (1:10)
                    Pines_write(2022, '7');
                end
                
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
     inst.SendScpi('SOUR:FUNC:MODE TASK');
     inst.SendScpi(':INIT:CONT ON');
     res = inst.SendScpi(':TRAC:DEL:ALL');
     assert(res.ErrCode==0);
end

function random_seq = get_random_seq(RMD_seq_length, seed)
    rng(seed);
    random_seq = randi([2, 3], 1, RMD_seq_length);
end

function rmd_seq = get_rmd_seq_block(pulse_t, n, t1, t2)
    %{
    Generate arbitrary length n-order Thue-Morse. t1 indicates + pulse. t2
    indicates b pulse.
    pulse_t = t1 or t2
    For example,
    get_rmd_seq_block(2,5,2,3) = \tilda{U_{5}}, get_rmd_seq_block(3,5,2,3) = U_{5}
    %}
    if n == 0
        rmd_seq = [pulse_t];
    elseif n > 0
        if pulse_t == t1
            rmd_seq = cat(2, get_rmd_seq_block(t1, n-1, t1, t2), get_rmd_seq_block(t2, n-1, t1, t2));
        elseif pulse_t == t2
            rmd_seq = cat(2, get_rmd_seq_block(t2, n-1, t1, t2), get_rmd_seq_block(t1, n-1, t1, t2));
        else
            error("pulse_t must be t1 or t2 (get_rmd_seq_block function)");
        end
    else
        error("input n must be integer or must be larger or equal to 0");
    end
end

function generateDTQCseqIQ(ch, amps, frequencies, lengths, phases, spacings, markers1, trigs, seg_idx_l, num_init_spin_lock)
%{
ch: channel RMD sequence is being applied
amps: amplitude of each pulses -- there are only three pulses pi/2
theta_a and theta_b
frequencies: always set to 0 for all pulses
lengths: length (in seconds) for each pulses
phases: phases of pulses in rotating frame
spacings: spacings after each pulse
markers1: always 1
trigs: triggers the digitizer
n_order: parameter for randomess n_order = 0: random, n_order = 1 dipole, n_order = 2 quadrupolar
%}
global inst
global sampleRateDAC
global interp
global sampleRateInterp 
    
    function setTask_Pulse(ch, seg_idx_l, num_init_spin_lock)
    %{
    4 types of segment:
    1: pi/2, 2: pi/2(1-0.1), 3: pi/2(1+0.1), 4: holding (64pts)
    %}
    disp('setting task table')
    % t_idx = tasktable index
    t_idx = 1;
    %% Initialize task table
    inst.SendScpi(sprintf(':INST:CHAN %d',ch));
    inst.SendScpi('TASK:ZERO:ALL');
    % set tasktable length to length of seg_idx_l + 4
    % +2 comes from the first and last holding segment and the +2 comes
    % from the first pi/2 y-pulse and num_init_spin_lock of pi/2 x-pulse
    inst.SendScpi(sprintf(':TASK:COMP:LENG %d', length(seg_idx_l)+4));
    %% set the first holding segment
    inst.SendScpi(':TASK:COMP:ENAB INT');
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
    t_idx = t_idx + 1;
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',5));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    %% set the initial pi/2 y-pulse
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
    t_idx = t_idx + 1;
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', 1));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', 1));
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    
    %% set the number pi/2 x-pulse to initialize state to +x
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
    t_idx = t_idx + 1;
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', 3));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', num_init_spin_lock));
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    
    %% set tasktable to apply random sequence
    for idx = 1: length(seg_idx_l)
        inst.SendScpi(sprintf(':TASK:COMP:SEL %d',t_idx));
        t_idx = t_idx + 1;
        inst.SendScpi(sprintf(':TASK:COMP:SEGM %d', seg_idx_l(idx, 1)));
        inst.SendScpi(sprintf(':TASK:COMP:LOOP %d', seg_idx_l(idx, 2)));
        inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',t_idx));
        inst.SendScpi(':TASK:COMP:TYPE SING');
    end
    
    %% set final segment to apply tasktable
    inst.SendScpi(sprintf(':TASK:COMP:SEL %d', t_idx));
    inst.SendScpi(sprintf(':TASK:COMP:LOOP %d',1));
    inst.SendScpi(sprintf(':TASK:COMP:SEGM %d',5));
    inst.SendScpi(':TASK:COMP:TYPE SING');
    inst.SendScpi(sprintf(':TASK:COMP:NEXT1 %d',1));
    inst.SendScpi('TASK:COMP:WRITE');
    resp = inst.SendScpi('SOUR:FUNC:MODE TASK');
    assert(resp.ErrCode==0);
    end

    function [I, Q, M1, Tr1] = get_square_pulse(frequency, lengthsPt, spacingsPt, amp, phase, marker, trig)
        %{
        Creates a square pulse that has following form.
        (I)
        _________
          |     |
          |     |
        (amp)   | 
          |     |
          |     |____________
        <------> <----------->
        (length)  (spacing)  
        
        (Q)
        _________
          |     |
          |     |
        (amp)   | 
          |     |
          |     |____________
        <------> <----------->
        (length)  (spacing)
        - ratio of (I) and (Q) controls the phase of the pulse in rotating
        frame.
        - M1(marker1) is turned on when the pulse is low (during spacing)
        - Tr1(marker2) is turned on only when reading out during spacing.
        %}
        DClen = spacingsPt;
        %%% make DC segment %%%
        [tempDCi, tempDCq] = makeDC(DClen);
        DClenreal = length(tempDCi);
        markDC = uint8(zeros(DClenreal, 1));
        markDC2 = uint8(zeros(DClenreal,1));
        %%% make Pulse %%% 
        [tempI, tempQ] = makeSqPulse(frequency, lengthsPt, amp, phase, 0, sampleRateDAC);
        pulseLenReal = length(tempI);
        markIQ = uint8(zeros(pulseLenReal, 1) + marker);
        markIQ2 = uint8(zeros(pulseLenReal, 1) + trig);
        I = [tempI tempDCi]; % first row is In-phase of pulse
        Q = [tempQ tempDCq]; % second row is Quadrature of pulse
        M1 = [markIQ' markDC']; % third row is blanking signal (M1)
        Tr1 = [markIQ2' markDC2'];
    end
    
    function segMat = assign_segMat(x, segMat, I, Q, M1, Tr1)
        segMat{1,x} = I;
        segMat{2,x} = Q;
        segMat{3,x} = M1;
        segMat{4,x} = Tr1;
    end

    %%%% FUNCTION STARTS HERE %%%%
    % Download all segments. +1 is there to create holding DC segment that
    % will be repeated.
    numSegs = 1 + length(lengths);
    lengthsPts = lengths * sampleRateDAC;
    spacingsPts = spacings * sampleRateDAC;
    
    disp('generating RMD RF pulse sequence')
    global segMat
    segMat = cell(4, numSegs); % added a fourth row for Marker2 (trigs)
    
    %% Make all segments
    for x = (1: 4)
        fprintf("%d \n", x);
        [I, Q, M1, Tr1] = get_square_pulse(frequencies(x), lengthsPts(x), spacingsPts(x), amps(x), phases(x), markers1(x), trigs(x));
        segMat = assign_segMat(x, segMat, I, Q, M1, Tr1);
    end
    
   
    
    %% Lastly, MAKE HOLDING DC SEGMENT
     x = numSegs;
    DClen = 64;
    [holdI, holdQ] = makeDC(DClen);
    markHold = uint8(zeros(DClen, 1));
    segMat = assign_segMat(x, segMat, holdI, holdQ, markHold, markHold);
    
    %% DOWNLOAD ALL SEGMENTS
    fprintf('pulse sequence generated')
    for x = (1: numSegs)
        downLoadIQ(ch, x, segMat{1,x}, segMat{2,x}, inst);
        downLoad_mrkr(ch, x, segMat{3,x}, segMat{4,x}, inst);
    end
    % NEED TO CHANGE PLANS
    % I cant be downloading segments for all the remainder -- seems like
    % memory error keeps happening after I download ~ 1600 segments
    % planning to do the heavy loading in the task table
    % simply going to repeat DC segment multiple times
    %% TODO: generate tasktable that can map content from 
    setTask_Pulse(ch, seg_idx_l, num_init_spin_lock);
 
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
 
 