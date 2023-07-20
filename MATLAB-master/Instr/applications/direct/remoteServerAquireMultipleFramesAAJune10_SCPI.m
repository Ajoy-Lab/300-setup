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
fullScaleMilliVolts =1000;
trigSource = 1; % 1 = external-trigger
dacChanInd = 1; % 1 = chan 1 and 2 = chan 2
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

dll_path = 'C:\Windows\System32\TEPAdmin.dll';

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
    
    inst.SendScpi(':DIG:CHAN:RANG HIGH');
    
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
                
%                 rc = inst.SetAdcCaptureTrigSource(adcChanInd, trigSource);
%                 assert(rc == 0);
                
                inst.SendScpi(':DIG:TRIG:SOUR EXT');
                
                inst.SendScpi(':DIG:TRIG:TYPE GATE');
                inst.SendScpi(':DIG:TRIG:SLOP NEG');
                
%                 rc = inst.SetAdcExternTrigPattern(0);
%                 assert(rc==0)
%                 
%                 rc = inst.SetAdcExternTrigGateModeEn(on);
%                 assert(rc==0);
                
                %
                inst.SendScpi(sprintf(':DIG:TRIG:LEV1 %g', 0.3));
%                 rc = inst.SetAdcExternTrigThresh(0,0.3);
%                 assert(rc==0);
                
                % rc = inst.SetAdcExternTrigDelay(0.0000025);
%                 rc = inst.SetAdcExternTrigDelay(delay+0.000001+33.8e-6+33.8e-6 - 12e-6);
                inst.SendScpi(sprintf(':DIG:TRIG:DEL:EXT %g', 4.1e-6));

%                 rc = inst.SetAdcExternTrigDelay(4.1e-6);
%                 assert(rc==0);

                %inst.SendScpi(':DIG:INIT ON');
                
%                 rc = inst.SetAdcAcquisitionEn(on,off);
%                 assert(rc == 0);
                
                fprintf('Instr setup complete and ready to aquire\n');
                
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
                
                netArray = NET.createArray('System.UInt16', readLen*numberOfPulses); %total array -- all memory needed
                
%                 rc = inst.SetAdcAcquisitionEn(on,off);
%                 assert(rc == 0);

                % Setup frames layout    
                inst.SendScpi(sprintf(':DIG:ACQ:FRAM:DEF %d, %d',numberOfPulses*loops, readLen));

%                 rc = inst.SetAdcFramesLayout(numberOfPulses*loops, readLen); %set memory of the AWG
%                 assert(rc == 0);
                
                fprintf('Waiting... Listen for Shuttle\n');
                inst.SendScpi(':DIG:INIT ON');
%                 rc = inst.SetAdcCaptureEnable(on);
%                 assert(rc == 0);
                
            case 3 % Measure
                
                %READ data from Sage and save
%                                         add this pause for measuring T1:
%                                         pause(postime*1e-3);
                                       % Sage_write('3');
                                       
                                       %KAH 4/19/22 insert case 3
                                       % Wait to see if AWG is ready to start capturing (?)
%                                        status = inst.ReadAdcCaptureStatus();
%                                        for i = 1 : 250
%                                            if status ~= 0
%                                                break;
%                                            end
%                                            pause(0.01);
%                                            status = inst.ReadAdcCaptureStatus();
%                                        end
% %                                        
%                                        rc = inst.SetAdcCaptureEnable(1);
%                                        assert(rc == 0);
% 
%                                        rc = inst.ReadAdcCaptureStatus();
                                       %assert(rc == 0);
%                                         disp('Press a key !')    
%                                         pause;
                                       % Set this to capture all data
                                       res = inst.SendScpi('DIG:ACQ:CAPT:ALL');
                                       assert(res.ErrCode == 0);
                                       
                                       % Other code has it set the number
                                       % of frames to capture
                                       cmd = [':DIG:ACQuire:FRAM:CAPT 1,' num2str(numberOfPulses*loops)]; %hard-coded in firstIndex=1
                                       %setting full data as total # frames
                                       res = inst.SendScpi(cmd);
                                       assert(res.ErrCode == 0);
                                       
                                       % Clean memory - note on 5/9/22 this
                                       % was taking too long so commented
                                       % out and memory cleaned afterwards.
                                       %res = inst.SendScpi(':DIG:ACQ:zero:ALL');
                                       %assert(res.ErrCode == 0);d
                                     
                                       % Close instrument in case improperly closed from prior run
                                       res = inst.SendScpi(':DIG:INIT OFF');
                                       assert(res.ErrCode == 0);
                                       % Starts the digitizer to wait for
                                       % capture
                                       res = inst.SendScpi(':DIG:INIT ON');
                                       assert(res.ErrCode == 0);
%                                        % Tells proteus to capture
%                                        % immediately on trig
%                                        res = inst.SendScpi(':DIG:TRIG:IMM');
%                                        assert(res.ErrCode == 0);
                                       fprintf('Waiting for trigger \n');
                                       rTrg=0;
                                       fprintf(['pulses expected =', num2str([numberOfPulses*loops]), '\n']);
                                       while rTrg < numberOfPulses*loops %set to total # triggers
                                           %res =
                                           %inst.SendScpi(':DIG:TRIG:IMM');
                                           res = inst.SendScpi(":DIG:ACQuire:FRAM:STATus?");
                                           assert(res.ErrCode == 0);
                                           respFields= split(convertCharsToStrings(char(res.RespStr)),',');
                                           %fprintf(1, '\nADC Status  ''%s''\n', char(res.RespStr));
                                           paramlen = size(respFields);
                                           if paramlen(1) >=4
                                               rTrg = str2num(respFields(4)); %the 4th element of respFields is numberOfPulses
                                           end
                                       end
                                       fprintf(1, '\nADC Status  ''%s''\n', char(res.RespStr));
                                       % Turn proteus acquisition off
                                       % before sending data to computer
                                       res = inst.SendScpi(':DIG:INIT OFF');
                                       assert(res.ErrCode == 0);
                                       
                                       % Choose what to read 
                                       % (only the frame-data without the header in this example)
                                       res = inst.SendScpi(':DIG:DATA:TYPE FRAM');
                                       assert(res.ErrCode == 0);
                                       fprintf('Transfering acquired data to computer....\n')
                                       
                                       pulseAmp = [];
                                       relPhase = [];
                                       %TESTpulses = [];
                                        for n = 1:loops
                                            fprintf('Start Read %d .... ', n);
                                            firstIndex = ((n-1)*numberOfPulses)+1;
                                            tic
                                            % Select frames starting at
                                            % next index
                                            cmd = [':DIG:DATA:SEL ' num2str(firstIndex) ',' num2str(numberOfPulses)];
                                            res = inst.SendScpi(cmd);
                                            assert(res.ErrCode == 0);
                                            
                                            rc = inst.ReadMultipleAdcFrames(adcDualChanMode, firstIndex, numberOfPulses, netArray); %this is where the device reads
                                            assert(rc == 0);
                                            samples = double(netArray); %get the data 
                                            fprintf('Read %d Complete\n', n);
                                            toc
                                            tic
                                            %delete mem
                                            fprintf('Clear mem %d .... ', n);
                                            rc =inst.WipeMultipleAdcFrames(adcDualChanMode, firstIndex, numberOfPulses, 0);
                                            assert(rc == 0);
                                            cmd = [':DIG:ACQ:ZERO ' num2str(firstIndex) ',' num2str(numberOfPulses) ',' num2str(0)];
                                            res = inst.SendScpi(cmd);
                                            assert(res.ErrCode == 0);
                                            fprintf('Clear mem %d Complete\n', n);
                                            toc
                                            
                                            tic
                                            fprintf('Starting iteration %d data processing....', n);
                                            pulses = reshape(samples, [], numberOfPulses); % reshape samples into a more useful array (2 dimensions)
                                            
%Save first 1s chunk of raw data 
%                                             savesinglechunk=true;
%                                             if savesinglechunk
%                                                 if  n==1 %determines which chunk will be saved
%                                                     
%                                                     pulsechunk = int16(pulses);
%                                                     a = datestr(now,'yyyy-mm-dd-HHMMSS');
%                                                     fn = sprintf([a,'_Proteus_chunk', num2str(n)]);
%                                                     Save data
%                                                     fprintf('Writing data to X:.....\n');
%                                                     lol=pulsechunk(:,2);
%                                                     save(['X:\' fn],'lol');
%                                                     writematrix(pulses);
%                                                     
%                                                 end
%                                             end
                                            
                                            %Pulses is now in an array that
                                            %has readLen # of rows and numberOfPulses # of columns
                                            %We iterate over each column (pulse):
%                                             for pulse_index=1:size(pulses,2)
%                                                 single_pulse_data=transpose(pulses(:,pulse_index));
%                                                 %FFT while subtracting DC component
%                                                 % Should we truncate some
%                                                 % data points to get rid of
%                                                 % first spike/ringing?
%                                                 Y=fft(single_pulse_data-mean(single_pulse_data));
%                                                 L=length(Y);
%                                                 P1 = Y(1:L/2+1); %Why the plus 1
%                                                 P1(2:end-1) = 2*P1(2:end-1); % This doubles the amplitude of everything but the first and last points. Why?
%                                                 f = sampleRate*(0:(L/2))/L;
%                                                 %start_fig(102,[2 2]);clf;plot_preliminaries(f/1e6,abs(P1),1,'nomarker');set(gca,'xlim',[0, 40]);
%                                                 %plot_labels('Frequency [MHz]','FT Signal[au]');
%                                                 %Grab out 20MHz component:
%                                                 [v,b]=min(abs(f-20e6));
%                                                 Signal(pulse_index)=P1(b);
%                                                 %TESTpulses = [TESTpulses, single_pulse_data];
%                                             end
%                                             dt=100e-6; %Approx amount of time between pulse/acquire windows
%                                             time_axis=0:dt:dt*(size(pulses,2)-1);
% %                                              start_fig(102,[2 2]);clf;
% %                                              p1=plot_preliminaries(time_axis,abs(Signal),1);set(p1,'markersize',2);
% %                                              plot_labels('Time [s]','Signal[au]');
                                            clear samples;
                                            
                                            %Set up variables to store data
                                            %in for QuantumStream
                                            %pulseAmp = [];
                                            %relPhase = []; %may need to move these outside the loop...
                                            
                                            % Do we want to include this
                                            % zero padding before above
                                            % data processing?
                                            power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
                                            padded_len= 2^(power_of_2) ;%2^15;
                                            % We already essentially do the
                                            % next few lines
                                            dF = sampleRate/padded_len; %set the discretization of freq in terms of sampleRate
                                            f = -sampleRate/2:dF:sampleRate/2-dF; %sampleRate sets the 'bandwidth'
% 
                                            [v,b1]=min(abs(f-20.00706e6)); %picks out the 20MHz component
                                            [v,b2]=min(abs(f+20.00706e6));
                                            % This truncates data but we
                                            % don't use?
                                            eff_read=100:readLen-100;
                                            cyclesPoints = 50;
                                            
                                            
                                            for i = 1:numberOfPulses
                                                 pulse = pulses(:, i);
%                                                  if n == 1
%                                                     if i == 2
%                                                         figure(4);clf;
%                                                         plot(pulse);
%                                                         xlabel('Time [\mus]');
%                                                         ylabel('Signal[a.u.]');
%                                                     end
%                                                  end
                                                pulseMean = mean(pulse);
                                                pulseDC = pulse - pulseMean; % remove DC
                                                X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                                                linFFT = (abs(X)/readLen);
                                                amp=(linFFT(b1) + linFFT(b2));
                                                phase1=angle(X(b1));
                                                phase2=angle(X(b2));
                                                %phase=phase1-phase2; %why does this exist?
                                                pulseAmp(i+(numberOfPulses*(n-1))) = amp;
                                                relPhase(i+(numberOfPulses*(n-1))) = phase1;
                                                %pulseAmp = [pulseAmp, amp];
                                                %relPhase = [relPhase, phase1];
                                            end
                                            save(['X:\d'] ,'pulses');
                                            
                                            clear pulses;
                                            fprintf('Data processing iteration %d complete!\n', n);
                                            toc
                                        end
                                        
                                        ivec=1:length(pulseAmp);
                                        delay2 = 0.000003; % dead time the unknown one, this is actually rof3 -Ozgur aka this is probably something we have to change?
                                        time_cycle = ((pw)+delay2+(tacq+2+4+10)*1e-6)*3 %this may be set differently as we plan to trigger on falling edge, time to do a pulse-acquire
                                        % using pulsed_spinlock_altskip, we take two time cycles to collect a data point
                                        time_axis=time_cycle.*ivec;
                                        %drop first point
                                        time_axis(1)=[];pulseAmp(1)=[];relPhase(1)=[]; %why only drop the first point? why drop it at all?
                                       
% % %                 inst.SendScpi(':DIG:INIT ON');
% % %                 
% % %                 for n = 1:250
% % %                     resp = inst.SendScpi(':DIG:ACQ:FRAM:STAT?');
% % %                     resp = strtrim(pfunc.netStrToStr(resp.RespStr));
% % %                     items = split(resp, ',');
% % %                     items = str2double(items);
% % %                     if length(items) >= 3 && items(2) == 1
% % %                         break
% % %                     end
% % %                     if mod(n, 10) == 0                
% % %                         fprintf('%d. %s Time:\n', fix(n / 10), resp);
% % %         %                 toc                
% % %                     end
% % %                     pause(0.1); 
% % %                 end
% % %                 
% % % %                 rc = inst.SetAdcCaptureEnable(on);
% % % %                 assert(rc == 0);
% % %                 
% % % %                 rc = inst.ReadAdcCaptureStatus();
% % %                 
% % %                 %                 % Wait until the capture completes
% % %                 %                 status = inst.ReadAdcCompleteFramesCount();
% % %                 %                 while status ~= numberOfPulses_total
% % %                 %                     pause(0.01);
% % %                 %                     status = inst.ReadAdcCompleteFramesCount();
% % %                 %                 end
% % %                 
% % %                 %                 readSize = uint64(readLen);
% % %                 %                 readOffset = uint64(offLen);
% % %                 
% % %                 chanIndex = 0;
% % %                 pulseAmp = [];
% % %                 relPhase = [];
% % %                 
% % %                 power_of_2 = floor(log2(readLen)); % this is normally 15 for 32us captures
% % %                 padded_len= 2^(power_of_2) ;%2^15;
% % %                 dF = sampleRate/padded_len; %set the discretization of freq in terms of sampleRate
% % %                 f = -sampleRate/2:dF:sampleRate/2-dF; %sampleRate sets the 'bandwidth'
% % %                 
% % %                 [v,b1]=min(abs(f-20.00706e6)); %picks out the 20MHz component
% % %                 [v,b2]=min(abs(f+20.00706e6));
% % %                 
% % %                 eff_read=100:readLen-100;
% % %                 cyclesPoints = 50;
% % %                 fprintf('Shuttle complete\n')
% % %                 fprintf('Transfering aquired data to computer....\n')
% % %                 for n = 1:loops
% % %                     fprintf('Start Read %d .... ', n);
% % %                     firstIndex = ((n-1)*numberOfPulses)+1;
% % %                     tic
% % %                     rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, numberOfPulses, netArray); %this is where the device reads
% % %                     assert(rc == 0);
% % %                     samples = double(netArray); %get the data (1s chunk)
% % %                     fprintf('Read %d Complete\n', n);
% % %                     toc
% % %                     
% % %                     tic
% % %                     %delete mem
% % %                     fprintf('Clear mem %d .... ', n);
% % %                     inst.SendScpi(':DIG:ACQ:ZERO:ALL');
% % %                     %rc =inst.WipeMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, 0);
% % %                     %assert(rc == 0);
% % %                     fprintf('Clear mem %d Complete\n', n);
% % %                     toc
% % %                     
% % %                     tic
% % %                     fprintf('Starting iteration %d data processing....', n);
% % %                     pulses = reshape(samples, [], numberOfPulses); % reshape samples into a more useful array (2 dimensions)
% % %                     
% % %                     if savealldata
% % %                         pulsechunk = int16(pulses);
% % %                         a = datestr(now,'yyyy-mm-dd-HHMMSS');
% % %                         fn = sprintf([a,'_Proteus','_chunk', num2str(n)]);
% % %                         % Save data
% % %                         fprintf('Writing data to Z:.....\n');
% % %                         save(['Z:\' fn],'pulsechunk');
% % %                         %writematrix(pulses);
% % %                     end
% % %                     
% % %                     if savesinglechunk
% % %                         if  n==1 %determines which chunk will be saved
% % %                             pulsechunk = int16(pulses);
% % %                             a = datestr(now,'yyyy-mm-dd-HHMMSS');
% % %                             fn = sprintf([a,'_Proteus_chunk', num2str(n)]);
% % %                             % Save data
% % %                             fprintf('Writing data to Z:.....\n');
% % %                             save(['Z:\' fn],'pulsechunk');
% % %                             %writematrix(pulses);
% % %                         end
% % %                     end
% % %                     
% % %                     clear samples;
% % %                     for i = 1:numberOfPulses
% % %                         pulse = pulses(:, i);
% % %                         %                         pulse = pulse(1024:readLen);
% % %                         %                         readLen=length(pulse);
% % %                         if savemultwind %save multiple consecutive windows of raw data if true
% % %                             if n==2 & i<=8
% % %                                 pulsechunk = int16(pulse);
% % %                                 a = datestr(now,'yyyy-mm-dd-HHMMSS');
% % %                                 fn = sprintf([a,'_Proteus_chunk', num2str(n) '_' num2str(i)]);
% % %                                 % Save data
% % %                                 fprintf('Writing data to Z:.....\n');
% % %                                 save(['Z:\' fn],'pulsechunk');
% % %                             end
% % %                         end
% % %                         
% % % %                         if n == 1
% % % %                             if i == 1
% % % %                                 figure(2);clf;
% % % %                                 plot(pulse);
% % % %                                 figure(3);clf;
% % % %                                 bandwidth=1/2/1e-9;
% % % %                                 Ttot=readLen*1e-9;
% % % %                                 F=linspace(-bandwidth,bandwidth,readLen);
% % % %                                 plot(F,abs(fftshift(fft(pulse))));
% % % %                                 hold on;
% % % %                                 yline(2048);
% % % %                             end
% % % %                         end
% % %                         if n == 1
% % %                             if i == 2
% % %                                 figure(4);clf;
% % %                                 plot(pulse);
% % % %                                 figure(5);clf;
% % % %                                 plot(F,abs(fftshift(fft(pulse))));
% % % %                                 hold on;
% % % %                                 yline(2048);
% % %                             end
% % %                         end
% % % %                         if n == 1
% % % %                             if i == 3
% % % %                                 figure(6);clf;
% % % %                                 plot(pulse);
% % % %                                 figure(7);clf;
% % % %                                 plot(F,abs(fftshift(fft(pulse))));
% % % %                                 hold on;
% % % %                                 yline(2048);
% % % %                             end
% % % %                         end
% % % % %                         pulse(1:41800)=[];
% % %                         pulseMean = mean(pulse);
% % %                         pulseDC = pulse - pulseMean; % remove DC
% % %                         X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
% % %                         linFFT = (abs(X)/readLen);
% % %                         %                         [amp, loc] = max(linFFT);
% % %                         amp=(linFFT(b1) + linFFT(b2));
% % %                         phase1=angle(X(b1));
% % %                         phase2=angle(X(b2));
% % %                         phase=phase1-phase2;
% % %                         %             linFFT_vec(:,i)=linFFT;
% % %                         pulseAmp(i+(numberOfPulses*(n-1))) = amp;
% % %                         relPhase(i+(numberOfPulses*(n-1))) = phase1;
% % %                     end
% % % %                     if n == 1
% % % %                         figure(8);clf;
% % % %                         plot(pulse);
% % % %                         hold on;
% % % %                         yline(2048);
% % % %                     end
% % %                     clear pulses;
% % %                     fprintf('Data processing iteration %d complete!\n', n);
% % %                     toc
% % %                 end
% % %                 
% % %                 %ivec=1:numberOfPulses*loops;
% % %                 ivec=1:length(pulseAmp);
% % %                 delay2 = 0.000003; % dead time the unknown one, this is actually rof3 -Ozgur
% % %                 
% % %                 %time_cycle=pw+96+(tacq+2+4+2+delay2)*1e-6;
% % %                 time_cycle=pw+delay2+(tacq+2+4+2)*1e-6
% % %                                  %time_cycle=pw+extraDelay+(4+2+2+tacq+17)*1e-6;
% % %                 time_axis=time_cycle.*ivec;
% % %                 %drop first point
% % %                 time_axis(1)=[];pulseAmp(1)=[];relPhase(1)=[];
% % %                 try
% % %                     figure(1);clf;
% % %                     %                plot(time_axis,(pulseAmp),'b-');hold on;
% % %                      plot(time_axis,pulseAmp,'r-');
% % %                     %plot(1:1:length(pulseAmp),pulseAmp,'r-');
% % %                     set(gca,'ylim',[0 max(pulseAmp)*1.1]);
% % %                     
% % %                     figure(12);clf;
% % %                      plot(time_axis,relPhase);
% % %                     %plot(1:1:length(relPhase),relPhase);
% % %                     
% % %                     %                 figure(3);clf;
% % %                     %                 plot(pulse);
% % %                     %                 hold on;
% % %                     %                 yline(2048);
% % %                 catch
% % %                     disp('Plot error occured');
% % %                 end
% % %                 
% % %                 %fn=dataBytes; %filename
% % %                 a = datestr(now,'yyyy-mm-dd-HHMMSS');
% % %                 fn = sprintf([a,'_Proteus']);
% % %                 % Save data
% % %                 fprintf('Writing data to Z:.....\n');
% % %                 save(['Z:\' fn],'pulseAmp','time_axis','relPhase');
% % %                 fprintf('Save complete\n');
                
            case 4 % Cleanup, save and prepare for next experiment
                inst.SendScpi(':DIG:INIT OFF');
%                 rc = inst.SetAdcCaptureEnable(off);
                assert(rc == 0);
                % Free the memory space that was allocated for ADC capture
                % Delete all wfm memory
                inst.SendScpi(':DIG:ACQ:ZERO:ALL');
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
                sweep_sigma = cmdBytes(6);
                symm = cmdBytes(7);
                srs_freq = cmdBytes(8); % should be 3.4e9
                srs_amp = cmdBytes(9);
                
                res = inst.SendScpi(['INST:CHAN ' num2str(dacChanInd)]); % select channel 2
                assert(res.ErrCode == 0);
                
                bits = 8;
                
                rampTime = 1/sweep_freq;
                fCenter = awg_center_freq - srs_freq;
                fStart = fCenter - 0.5*awg_bw_freq;
                fStop = fCenter + 0.5*awg_bw_freq;
                dt = 1/sampleRateDAC;
                dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
                fprintf('waveform length - ');
                fprintf(num2str(length(dacSignal)));
                fprintf('\n') ;
                
                % Define segment 1 
                res = inst.SendScpi(strcat(':TRAC:DEF 1,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
                assert(res.ErrCode == 0);
                
                % select segmen 1 as the the programmable segment
                res = inst.SendScpi(':TRAC:SEL 1');
                assert(res.ErrCode == 0); 
                
                % Download the binary data to segment 1
                %res = inst.WriteBinaryData(':TRAC:DATA 0,#', dacSignal);
                dacSignal = uint8(dacSignal);
                res = inst.WriteBinaryData(':TRAC:DATA 0,', dacSignal);
                assert(res.ErrCode == 0);
                
                srs_freq_str = [':SOUR:NCO:CFR1 ' sprintf('%0.2e', srs_freq)];
                res = inst.SendScpi(srs_freq_str);
                assert(res.ErrCode == 0);
                
                inst.SendScpi(':NCO:SIXD1 ON');
                
                inst.SendScpi(':SOUR:MODE NCO');
                                
%                 res = inst.SendScpi(':SOUR:MODE IQM'); % IQ MODULATOR --
% %                 THINK OF A MIXER, BUT DIGITAL
%                 assert(res.ErrCode == 0);
                
                try
                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);
                catch
                sampleRateDAC_str = [':FREQ:RAST ' sprintf('%0.2e', sampleRateDAC)];
                res = inst.SendScpi(sampleRateDAC_str); % set sample clock
                assert(res.ErrCode == 0);
                end
                
                % Define segment 2
                res = inst.SendScpi(strcat(':TRAC:DEF 2,',num2str(length(dacSignal)))); % define memory location 1 of length dacSignal
                assert(res.ErrCode == 0);
                
                % select segmen 2 as the the programmable segment
                res = inst.SendScpi(':TRAC:SEL 2');
                assert(res.ErrCode == 0); 
                
                % Download the binary data to segment 2
                res = inst.WriteBinaryData(':TRAC:DATA 0,#', fliplr(dacSignal));
                assert(res.ErrCode == 0);
                
                srs_freq_str = [':SOUR:NCO:CFR1 ' sprintf('%0.2e', srs_freq)];
                res = inst.SendScpi(srs_freq_str);
                assert(res.ErrCode == 0);
                
                inst.SendScpi(':NCO:SIXD1 ON');
                
                inst.SendScpi(':SOUR:MODE NCO');
                
%                 res = inst.SendScpi(':SOUR:MODE IQM'); % IQ MODULATOR --
% %                 THINK OF A MIXER, BUT DIGITAL
%                 assert(res.ErrCode == 0);
                
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
                
                segLen = 40960;
                segLenDC = 50048;
                cycles = 800;
                
%                 dacSignalMod = sine(cycles, 0, segLen, bits); % sclk, cycles, phase, segLen, amp, bits, onTime(%)
%                 
%                 dacSignalDC = [];
%                 for i = 1 : segLenDC
%                   dacSignalDC(i) = 127;
%                 end
%                 
%                 dacSignal = [dacSignalMod dacSignalDC];
                
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
                
%                 amp_str = [':SOUR:VOLT ' sprintf('%0.2f', awg_amp)];
%                 res = inst.SendScpi(amp_str);
%                 assert(res.ErrCode == 0);
                
%                 res = inst.SendScpi(':OUTP ON');
%                 assert(res.ErrCode == 0);
                
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

fclose(u2);

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
