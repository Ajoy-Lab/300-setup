%% Clear everything
clear;
close;

%% Set defaults Vars
savealldata=false;
savesinglechunk=false;
savemultwind=false;

sampleRate = 1000e6;
sampleRateDAC = 9e9;
adcDualChanMode = 1;
% fullScaleMilliVolts =1000;
trigSource = 1; % 1 = external-trigger
dacChanInd = 2; % 1 = chan 1 and 2 = chan 2
adcChanInd = 2; % ADC Channel 1
measurementTimeSeconds = 7; %Integer
%delay = 0.0000178; % dead time
%delay = 0.0000348; % dead time
%delay=0.00000148;
%delay = 0.0000108; % dead time
%delay = 0.0000028; % dead time
delay = 0.0000038; % dead time


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
    
    % ---------------------------------------------------------------------
    % ADC Config
    % ---------------------------------------------------------------------
    
    inst.SendScpi(':DIG:MODE DUAL');
    
    inst.SendScpi(sprintf(':DIG:CHAN %d', adcChanInd)); 
   
    inst.SendScpi(sprintf(':DIG:FREQ %g', sampleRate));
    
    inst.SendScpi(':DIG:CHAN:RANG LOW');
    
    % Enable acquisition in the digitizer's channels  
    inst.SendScpi(':DIG:CHAN:STATE ENAB');
    
    fprintf('ADC Configured\n');
    
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
    
    % Loop to repeatedly wait for messages and send replies
    % Break or Ctrl+C to get out of loop
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
                

                
            case 2 % Aquire on trig
                
                fprintf('Calculate and set data structures...\n');
                
                   pw = cmdBytes(2)*1e-6;
%                 pw=110*1e-6;
%                 pw= 30.0*1e-6
                %pw=32.5e-6
               pw=45.0e-6;
%                  
                   %pw=85e-6;
%                   pw=68e-6;
%                   pw=42.5e-6;
%                   pw=34e-6;
%                   pw=136e-6;
%                   pw=170e-6;
%                   %pw=113.333e-6;
%                   %pw=56.666e-6;
                %pw=40.0e-6;
                
                numberOfPulses_total= cmdBytes(3);
%                 numberOfPulses_total= 113636;
%                 numberOfPulses_total= 12000;
%                     numberOfPulses_total= 234375;
%                 numberOfPulses_total= cmdBytes(3)*floor(8*1.5^cmdBytes(2))+6000;
%                 numberOfPulses_total=169993;
                
                Tmax=cmdBytes(4);
                
                
                tacq=cmdBytes(5);
%                 tacq=128;
%                 tacq=64;
%                 tacq=96;
                
                numberOfPulses= floor(numberOfPulses_total/Tmax); %in 1 second %will be 1 for FID
                loops=Tmax;
                
                 readLen = round2((tacq+2)*1e-6/1e-9,96)-96;
                 %readLen = round2((tacq)*1e-6/1e-9,96)-96;
%                 readLen = 33888; % round2((tacq+2)*1e-6/1e-9,96)-96 %constraint: has to be multiple of 96, add 4 of deadtime %number of points in a frame
                  %readLen = 16896; % for tacq=16
%                 readLen = 65856; % for tacq=64
%                 readLen = 129888; % for tacq=128

%                 readLen = 66912; % actual tacq=64
%                 readLen = 131040;
%                 readLen = 129600;
%                 readLen = 126912; % actual tacq=128
%                 readLen = 97920; % for tacq=96
                
                offLen = 0;
                
                %                 rc = inst.SetAdcCaptureTrigSource(adcChanInd, trigSource);
%                 assert(rc == 0);
                
                 rc = inst.SendScpi(sprintf(':DIG:ACQ:FRAM:DEF %d, %d',numberOfPulses*loops, readLen));
                 assert(rc.ErrCode == 0);
%                 inst.SendScpi(':DIG:TRIG:SOUR EXT');
%                 
%                 inst.SendScpi(':DIG:TRIG:TYPE GATE');
%                 inst.SendScpi(':DIG:TRIG:SLOP NEG');
                
                 rc = inst.SendScpi(':DIG:TRIG:SOUR EXT'); %digChan
                 assert(rc.ErrCode == 0);
                 rc = inst.SendScpi(':DIG:TRIG:SLOP NEG');
                 assert(rc.ErrCode == 0);
                %inst.SendScpi(sprintf(':DIG:TRIG:SELF %f', 0.0)); %0.025 
                
%                 rc = inst.SetAdcExternTrigPattern(0);
%                 assert(rc==0)
%                 
%                 rc = inst.SetAdcExternTrigGateModeEn(on);
%                 assert(rc==0);
                
                %
                rc = inst.SendScpi(sprintf(':DIG:TRIG:LEV1 %g', 1));
                assert(rc.ErrCode == 0);
%                 rc = inst.SetAdcExternTrigThresh(0,0.3);
%                 assert(rc==0);
                
                % rc = inst.SetAdcExternTrigDelay(0.0000025);
%                 rc = inst.SetAdcExternTrigDelay(delay+0.000001+33.8e-6+33.8e-6 - 12e-6);
                rc = inst.SendScpi(sprintf(':DIG:TRIG:DEL:EXT %g', 4.1e-6)); 
                assert(rc.ErrCode == 0);

%                 rc = inst.SetAdcExternTrigDelay(4.1e-6);
%                 assert(rc==0);

                %inst.SendScpi(':DIG:INIT ON');
                
%                 rc = inst.SetAdcAcquisitionEn(on,off);
%                 assert(rc == 0);
                
                fprintf('Instr setup complete and ready to aquire\n');
                
                netArray = NET.createArray('System.UInt16', readLen*numberOfPulses); %total array -- all memory needed
                
%                 rc = inst.SetAdcAcquisitionEn(on,off);
%                 assert(rc == 0);

                % Setup frames layout    

                
                 rc = inst.SendScpi(':DIG:ACQ:FRAM:CAPT:ALL');   
                 assert(rc.ErrCode == 0);
                 %rc = inst.SendScpi(':DIG:ACQ:ZERO:ALL');
                 assert(rc.ErrCode == 0);
    %                 rc = inst.SetAdcFramesLayout(numberOfPulses*loops, readLen); %set memory of the AWG
    %                 assert(rc == 0);
%                  rc = resp = inst.SendScpi(':DIG:DATA:FORM?');
%                  assert(rc == 0);
%                  resp = strtrim(pfunc.netStrToStr(resp.RespStr));
    
                fprintf('Waiting... Listen for Shuttle\n');
                rc = inst.SendScpi(':DIG:INIT OFF'); 
                assert(rc.ErrCode == 0);
                rc = inst.SendScpi(':DIG:INIT ON');
                assert(rc.ErrCode == 0);
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
                
                netArray = NET.createArray('System.UInt16', readLen*numberOfPulses); %total array -- all memory needed
                                
                power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
                padded_len= 2^(power_of_2) ;%2^15;
                dF = sampleRate/padded_len; %set the discretization of freq in terms of sampleRate
                f = -sampleRate/2:dF:sampleRate/2-dF; %sampleRate sets the 'bandwidth'
                
                [v,b1]=min(abs(f-20.00706e6)); %picks out the 20MHz component
                [v,b2]=min(abs(f+20.00706e6));
                
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
                    
%                     rc = GetDigitizerData2(   inst,...
%                                         pfunc,...
%                                         cType,...
%                                         true,...
%                                         chanIndex,...
%                                         3,...
%                                         false,...                                   
%                                         0,...
%                                         0.0,...
%                                         sampleRate,...
%                                         0.0,...
%                                         1,...
%                                         readLen,...
%                                         numberOfPulses,...
%                                         0);
                                                       
                    rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, numberOfPulses, netArray); %this is where the device reads
%                     rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, nFrames, netArray2); %this is where the device reads
                    
                    assert(rc == 0);
%                     samples = double(netArray2); %get the data (1s chunk)
                    samples = double(netArray); %get the data (1s chunk)
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
                        
%                         if n == 1
%                             if i == 1
%                                 figure(2);clf;
%                                 plot(pulse);
%                                 figure(3);clf;
%                                 bandwidth=1/2/1e-9;
%                                 Ttot=readLen*1e-9;
%                                 F=linspace(-bandwidth,bandwidth,readLen);
%                                 plot(F,abs(fftshift(fft(pulse))));
%                                 hold on;
%                                 yline(2048);
%                             end
%                         end
                        if n == 1
                            if i == 2
                                figure(4);clf;
                                plot(pulse);
%                                 figure(5);clf;
%                                 plot(F,abs(fftshift(fft(pulse))));
%                                 hold on;
%                                 yline(2048);
                            end
                        end
%                         if n == 1
%                             if i == 3
%                                 figure(6);clf;
%                                 plot(pulse);
%                                 figure(7);clf;
%                                 plot(F,abs(fftshift(fft(pulse))));
%                                 hold on;
%                                 yline(2048);
%                             end
%                         end
% %                         pulse(1:41800)=[];
                        pulseMean = mean(pulse);
                        pulseDC = pulse - pulseMean; % remove DC
                        X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                        linFFT = (abs(X)/readLen);
                        %                         [amp, loc] = max(linFFT);
                        amp=(linFFT(b1) + linFFT(b2));
                        phase1=angle(X(b1));
                        phase2=angle(X(b2));
                        phase=phase1-phase2;
                        %             linFFT_vec(:,i)=linFFT;
                        pulseAmp(i+(numberOfPulses*(n-1))) = amp;
                        relPhase(i+(numberOfPulses*(n-1))) = phase2; %positive phase
                        %relPhase(i+(numberOfPulses*(n-1))) = phase1; %negative phase
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
                
                %ivec=1:numberOfPulses*loops;
                ivec=1:length(pulseAmp);
                delay2 = 0.000003; % dead time the unknown one, this is actually rof3 -Ozgur
                
                %time_cycle=pw+96+(tacq+2+4+2+delay2)*1e-6;
                time_cycle=pw+delay2+(tacq+2+4+2)*1e-6
                                 %time_cycle=pw+extraDelay+(4+2+2+tacq+17)*1e-6;
                time_axis=time_cycle.*ivec;
                %drop first point
                time_axis(1)=[];pulseAmp(1)=[];relPhase(1)=[];
                try
                    figure(1);clf;
                    %                plot(time_axis,(pulseAmp),'b-');hold on;
                     plot(time_axis,pulseAmp,'r-');
                    %plot(1:1:length(pulseAmp),pulseAmp,'r-');
                    set(gca,'ylim',[0 max(pulseAmp)*1.1]);
                    
                    figure(12);clf;
                     plot(time_axis,relPhase);
                    %plot(1:1:length(relPhase),relPhase);
                    
                    %                 figure(3);clf;
                    %                 plot(pulse);
                    %                 hold on;
                    %                 yline(2048);
                catch
                    disp('Plot error occured');
                end
                
                %fn=dataBytes; %filename
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                % Save data
                fprintf('Writing data to Z:.....\n');
                save(['Z:\' fn],'pulseAmp','time_axis','relPhase');
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
                
                
                
                
                
                
                awg_center_freq = cmdBytes(2);
                awg_bw_freq = cmdBytes(3);
                awg_amp = cmdBytes(4);
                sweep_freq = cmdBytes(5);
%                 sweep_freq = 3e3;
                sweep_sigma = cmdBytes(6);
                symm = cmdBytes(7);
                srs_freq = cmdBytes(8); % should be 3.4e9
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
                
                inst.SendScpi("*CLS; *RST")
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                bits = 8;
                
                rampTime = 1/sweep_freq;
                fCenter = awg_center_freq - srs_freq;
                fStart = fCenter - 0.5*awg_bw_freq;
                disp(['fstart = ' num2str(fStart)]);
                fStop = fCenter + 0.5*awg_bw_freq;
                dt = 1/sampleRateDAC;
                
%                 dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits); 
                
                chirps{1}.dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
                chirps{2}.dacSignal = fliplr(chirps{1}.dacSignal);
                chirps{1}.segLen = length(chirps{1}.dacSignal);
                chirps{2}.segLen = length(chirps{2}.dacSignal);
                chirps{1}.segm = 1;
                chirps{2}.segm = 2;
                chirps{1}.srs_freq = 0.3625e9;
                chirps{2}.srs_freq = 0.3625e9;
                
                fprintf('waveform length - ');
%                 fprintf(num2str(length(dacSignal)));
                fprintf('\n') ;
                
                
                
                pol_times = [];
                task_list = build_tasktable(inst,pol_times,chirps,sampleRateDAC);
                
                segLen = 40960;
                segLenDC = 50048;
                cycles = 800;

                
            case 7 % Play MW chirp waveform
                
                % ---------------------------------------------------------------------
                % Play segment 1 in channel 1
                % ---------------------------------------------------------------------
                
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':SOUR:FUNC:MODE:SEG 1');
                assert(res.ErrCode == 0);
                
                amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
                res = inst.SendScpi(amp_str);
                assert(res.ErrCode == 0);
                
                res = inst.SendScpi(':OUTP ON');
                assert(res.ErrCode == 0);
                
                fprintf('Waveform generated and playing\n');
                                fprintf('still playing\n');
                
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
                
                amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
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

function rNb = roundCheck(nb, units)
    global debug
    rNb = round(nb/units);
    if ~isempty(debug) && debug > 1 && abs(rNb*units - nb) > units/100
       disp([inputname(1) ' = ' num2str(nb/units) ' rounded to ' num2str(rNb)]);
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

        % Download the binary data to segment
        res = inst.WriteBinaryData(':TRAC:DATA 0,', chirps{i}.dacSignal);
        assert(res.ErrCode == 0);
        
        srs_freq_str = [':SOUR:NCO:CFR1 ' sprintf('%0.2e',  sampleRateDAC - chirps{i}.srs_freq)]; %srs_freq
        res = inst.SendScpi(srs_freq_str);
        assert(res.ErrCode == 0);

       res = inst.SendScpi(':NCO:SIXD1 ON');
       assert(res.ErrCode == 0);

        rc = inst.SendScpi(':SOUR:MODE DUC');
        assert(rc.ErrCode == 0);

        try
            sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
            res = inst.SendScpi(sampleRateDAC_str); % set sample clock
            assert(res.ErrCode == 0);
        catch
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
    task_list = {};
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
