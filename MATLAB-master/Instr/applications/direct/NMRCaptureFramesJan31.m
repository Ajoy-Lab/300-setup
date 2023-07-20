%% NMR Spin-Lock Experiment
% Capture Data with ADC

clear all;
close all;
clc;

%% Set defaults Vars
sampleRateADC = 1e9;
adcDualChanMode = 1;
fullScaleMilliVolts = 1000;
trigSource = 1; % 1 = external-trigger
chanIndex = 0; % ADC Channel 1

delay = 0.0000172; % minimum delay between pulse and acquisition
framesInOneWindow = 12499;

remoteAddr = '192.168.1.2';
remotePort = 2020;
localPort = 9090;

off = 0;
on = 1;

%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Windows\System32\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*

admin = CProteusAdmin(@OnLoggerEvent);
rc = admin.Open();
assert(rc == 0);

%% Use the administrator, and close it at the end
try
    slotIds = admin.GetSlotIds();
    numSlots = length(slotIds);
    assert(numSlots > 0);
    
    connect = 1;
    % If there are multiple slots, let the user select one ..
    sId = slotIds(1);
    if numSlots > 1
        fprintf(1, '\n%d slots were found\n', numSlots);
        for n = 1:numSlots
            sId = slotIds(n);
            slotInfo = admin.GetSlotInfo(sId);
            if ~slotInfo.IsSlotInUse
                modelName = slotInfo.ModelName;
                if slotInfo.IsDummySlot
                    fprintf(1, ' * Slot Number: Model %s [Dummy Slot].\n', sId, ModelName);
                else
                    fprintf(1, ' * Slot Number: Model %s.\n', sId, ModelName);
                end
            end
        end
        choice = input('Enter SlotId ');
        fprintf(1, '\n');
        sId = uint32(choice);
    end
    
    % Connect to the selected instrument ..
    should_reset = true;
    inst = admin.OpenInstrument(sId);%, should_reset);
    instId = inst.InstrId;
    %% ---------------------------------------------------------------------
    %  Setup instrument
    %  ---------------------------------------------------------------------
    
    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
    
    pause(0.01);
    
    res = inst.SendScpi('*CLS'); % clear
    assert(res.ErrCode == 0);
    
    res = inst.SendScpi('*RST'); % reset
    assert(res.ErrCode == 0);
    
%     res = inst.SendScpi(':FREQ:RAST 1E9'); % set sample clock
%     assert(res.ErrCode == 0);    
    
    fprintf('Reset complete\n');
    
    % ---------------------------------------------------------------------
    % Operate the ADC with direct functions 
    % ---------------------------------------------------------------------
    
    % Turn on ADC dual-channels mode (state = 1)
    rc = inst.SetAdcDualChanMode(on);
    assert(rc == 0);

    rc = inst.SetAdcSamplingRate(sampleRateADC);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(chanIndex,fullScaleMilliVolts); % mv - 500, 800, 1000
    assert(rc == 0);
    
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
%         Wait for message
        while(u2.BytesAvailable == 0)
%             If no bytes in u2 buffer, wait 10ms then check again
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
%         dataBytes=strtok(dataBytes,',');
        cmdByte=cmdBytes(1);
    
        % 1 - Set triggers
        % 2 - Initialize data structures
        % 3 - Acquire, read, process and write data to NAS
        % 4 - Cleanup and prepare for next cycle
        % 5 - Shut down and disconnect intrument
        
        switch cmdByte
            case 2 % Initialize data structures
                
                pw = cmdBytes(2)*1e-6;
                totalNumberOfFrames = cmdBytes(3);
                Tmax = cmdBytes(4);
                tacq = cmdBytes(5);
                
                readLen = round2((tacq+2)*sampleRateADC*1e-6,96)-96; %constraint: has to be multiple of 96, add 4 of deadtime
                numberOfWindows = floor(totalNumberOfFrames/framesInOneWindow);
                netArray = NET.createArray('System.UInt16', readLen*framesInOneWindow);
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);
                
                rc = inst.SetAdcFramesLayout(totalNumberOfFrames, readLen);
                assert(rc == 0);
                
                fprintf('ReadLen: %d ... Number of frames: %d\n', readLen, totalNumberOfFrames);
                fprintf('Capture size: %d \n', totalNumberOfFrames*readLen);
                
                rc = inst.SetAdcCaptureEnable(on);
                assert(rc == 0);
                
                fprintf('Aquisition enabled, starting measurement...\n');
                
            case 1 % Set triggers
                
%                 fprintf('Trigger set - ');
%                 trigSource = 1; % 2 = Self, 1 = External
                trigSource = 1;%input('Select trigger source - enter 1 for External or 2 for Self :');
                
                rc = inst.SetAdcCaptureTrigSource(chanIndex, trigSource); % 1 means external trigger, 2 means self trigger
                assert(rc == 0);
                
                if trigSource == 1
                    
                    rc = inst.SetAdcExternTrigPattern(0);
                    assert(rc==0);
                    
                    %Enable trigger gate mode
                    rc = inst.SetAdcExternTrigGateModeEn(1);
                    assert(rc==0);
                    
                    %Set trigger polarity 0-positive 1-negative
                    rc = inst.SetAdcExternTrigPolarity(1);
                    assert(rc==0);
                    
                    %Set trigger level to 1V
                    rc = inst.SetAdcExternTrigThresh(0,1);
                    assert(rc==0);
                    
                    %Set trigger delay
                    rc = inst.SetAdcExternTrigDelay(delay);
                    assert(rc==0);
                    
                end
                
                fprintf('Instrument setup complete and ready to aquire!\n');
                
            case 3 % Acquire, read, process and write data to NAS
                
                % Wait until the capture completes
                status = inst.ReadAdcCompleteFramesCount();
                
                while status ~= totalNumberOfFrames
                    pause(0.01);
                    status = inst.ReadAdcCompleteFramesCount();
                end
%                 
%                 prompt = 'Wait for shuttle to complete ';
%                 x = input(prompt);    
%                 fprintf('Reading Frames....\n'); 
                                
                totalAmp = [];
                power_of_2 = ceil(log2(readLen)); % this is normally 15 for 32us captures
                padded_len= 2^(power_of_2) ;%2^15;
                dF = sampleRateADC/padded_len; %set the discretization of freq in terms of sampleRate
                f = -sampleRateADC/2:dF:sampleRateADC/2-dF; %sampleRate sets the 'bandwidth'
                
                [v,b1]=min(abs(f-20e6)); %picks out the 20MHz component
                [v,b2]=min(abs(f+20e6));
                
                for n = 1 : numberOfWindows
                    
%                     Read one window of frames simultaneously
                    tic
                    fprintf('Start Read Window %d ....', n);
                    firstIndex = ((n-1)*framesInOneWindow)+1;
                    rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, framesInOneWindow, netArray); %this is where the device reads
                    assert(rc == 0);
                    samples = double(netArray); %get the data (1s chunk)
                    fprintf('Read Window %d Complete\n', n);
                    toc
                    
%                     Clear netArray by filling it back up with zeroes
                    tic
                    fprintf('Clear mem %d ....', n);
                    rc =inst.WipeMultipleAdcFrames(chanIndex, firstIndex, framesInOneWindow, 0);
                    assert(rc == 0);
                    fprintf('Clear mem %d Complete\n', n);
                    toc
                    
%                     Process data window one frame at a time
                    tic
                    fprintf('Starting data processing window %d ... ', n);
                    pulseWindow = reshape(samples, [], framesInOneWindow);
                    clear samples;
                    for i = 1 : framesInOneWindow
                        pulse = pulseWindow(:,i);
                        pulseMean = mean(pulse);
                        pulseDC = pulse - pulseMean; % remove DC
                        X = fftshift(fft(pulseDC,padded_len)); % perform FFT and zero-pad to the padded_len
                        linFFT = (abs(X)/readLen);
                        amp=(linFFT(b1) + linFFT(b2));
                        totalAmp(i+((n-1)*framesInOneWindow)) = amp; % save amplitude as one point
                    end
                    clear pulseWindow;
                    fprintf('Data processing window %d complete\n', n);
                    toc
                    
%                     Plot the last pulse of each window
                    figure(2);clf;
                    plot(pulseDC);
                    hold on;
                end
                
                ivec=1:framesInOneWindow*numberOfWindows;
                time_cycle=pw*1e-6+(4+2+2+tacq)*1e-6; % fix this ! !
                time_axis=time_cycle.*ivec;
                time_axis(1)=[];totalAmp(1)=[]; % drop first point (false trigger)
                
                try
                figure(1);clf;
                plot(time_axis,(smooth(totalAmp,10)),'r-');
                set(gca,'ylim',[0 max(totalAmp)*1.1]);
                catch
                    disp('Plot error occurred');
                end
                
%                 Save data
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                fprintf('Writing data to NAS ...\n');
                save(['Z:\' fn],'pulseAmp','time_axis');
                
            case 4 % Clean up and prepare for next experiment
                
                rc = inst.SetAdcCaptureEnable(off);
                assert(rc == 0);
                
                % Free the memory space that was allocated for ADC capture
                rc = inst.FreeAdcReservedSpace();
                assert(rc == 0);
                
                fprintf('Complete and ready for next acquisition!\n');
                
            case 5 % Shut down and disconnect instrument
                
                connect = 0;
                
        end % kill this switch
        
    end % end while
    
% ---------------------------------------------------------------------
% End of experiment
% ---------------------------------------------------------------------
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd - server stopped! \nInstrument Error Status = %s\n', netStrToStr(res.RespStr));
    
    % It is recommended to disconnect from instrumet at the end
    rc = admin.CloseInstrument(instId);
    
    % Close the administrator at the end ...
    admin.Close();
catch ME
    admin.Close();
    rethrow(ME)
end

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
