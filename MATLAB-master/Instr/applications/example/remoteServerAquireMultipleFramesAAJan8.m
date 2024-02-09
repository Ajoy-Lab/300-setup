%% Clear everything
clear;
close;

%% Set defaults Vars
sampleRate = 1000e6;
adcDualChanMode = 1;
fullScaleMilliVolts =1000;
trigSource = 1; % 1 = external-trigger
adcChanInd = 0; % ADC Channel 1
measurementTimeSeconds = 7; %Integer

remoteAddr = '192.168.1.2';
remotePort = 2020;
localPort = 9090;

off = 0;
on = 1;


%% Load TEPAdmin.dll which is a .Net Assembly
% The TEPAdmin.dll is installed by WDS Setup in C:\Windows\System32 

% Currently the tested DLL is installed in the following path
asm = NET.addAssembly('C:\Users\Simon\Documents\Tabor Electronics\ProteusAwg-10243-Customer\TEPAdmin.dll');

import TaborElec.Proteus.CLI.*
import TaborElec.Proteus.CLI.Admin.*
%% Create Administrator
% Create instance of the CProteusAdmin class, and open it.
% Note that only a single CProteusAdmin instance can be open, 
% and it must be kept open till the end of the session.

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
    should_reset = false;
    inst = admin.OpenInstrument(sId);
    instId = inst.InstrId;

    % ---------------------------------------------------------------------
    % Setup instrument
    % ---------------------------------------------------------------------

    res = inst.SendScpi('*IDN?');
    assert(res.ErrCode == 0);
    fprintf(1, '\nConnected to ''%s''\n', netStrToStr(res.RespStr));
   
    pause(0.01);
    
    
    res = inst.SendScpi('*CLS');
    assert(res.ErrCode == 0);

    res = inst.SendScpi('*RST');
    assert(res.ErrCode == 0);

    fprintf('Reset complete\n');

    % ---------------------------------------------------------------------
    % ADC Config
    % ---------------------------------------------------------------------
    
    rc = inst.SetAdcDualChanMode(1); %set to dual frame granularity = 48
    assert(rc == 0);
    
    rc = inst.SetAdcSamplingRate(sampleRate);
    assert(rc == 0);
    
    rc = inst.SetAdcFullScaleMilliVolts(0,fullScaleMilliVolts);
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
        %cmdByte
        % 1 - Init
        % 2 - Aquire on Trig
        % 3 - Measure
        % 4 - Cleanup - prepare for next cycle
        % 5 - Disconnect Intrument
        
        %measure_time=str2double(dataBytes); %in sec
        switch cmdByte
            case 1 % Setup
                
                rc = inst.SetAdcCaptureTrigSource(adcChanInd, trigSource);
                assert(rc == 0);
                
                rc = inst.SetAdcExternTrigPattern(0);
                assert(rc==0)
                
                rc = inst.SetAdcExternTrigGateModeEn(on);
                assert(rc==0);
                
                rc = inst.SetAdcExternTrigThresh(0,0.3);
                assert(rc==0);
                
                %rc = inst.SetAdcExternTrigDelay(0.0000025);
                rc = inst.SetAdcExternTrigDelay(0.0000967);
                assert(rc==0); 
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);

                fprintf('Instr setup complete and ready to aquire\n');
                
            case 2 % Aquire
 
                fprintf('Calculate and set data structures...\n');

                pw= cmdBytes(2)*1e-6;
                numberOfPulses_total= cmdBytes(3);
                Tmax=cmdBytes(4);
                tacq=cmdBytes(5);
                numberOfPulses= floor(numberOfPulses_total/Tmax); %in one second
                loops=Tmax;
                
                readLen = round2((tacq+2)*1e-6/1e-9,96)-96; %constaint: has to be multiple of 96, add 4 of deadtime
                
                offLen = 0;

                netArray = NET.createArray('System.UInt16', readLen*numberOfPulses);
                
                rc = inst.SetAdcAcquisitionEn(on,off);
                assert(rc == 0);
                
                rc = inst.SetAdcFramesLayout(numberOfPulses*loops, readLen);
                assert(rc == 0);

                fprintf('Waiting... Listen for Shuttle\n');
                rc = inst.SetAdcCaptureEnable(on);
                assert(rc == 0);
  
            case 3
                
                status = inst.ReadAdcCaptureStatus();
                for i = 1 : 250
                    if status ~= 0
                        break;
                    end
                    pause(0.01);
                    status = inst.ReadAdcCaptureStatus();
                end
                
                rc = inst.SetAdcCaptureEnable(on);
                assert(rc == 0);

                rc = inst.ReadAdcCaptureStatus();

                readSize = uint64(readLen);
                readOffset = uint64(offLen);

                chanIndex = 0;
                pulseAmp = [];
                padded_len=2^15;
                dF = sampleRate/padded_len;
                f = -sampleRate/2:dF:sampleRate/2-dF;

                [v,b1]=min(abs(f-20e6));
                [v,b2]=min(abs(f+20e6));

                eff_read=100:readLen-100;
                cyclesPoints = 50;
                fprintf('Shuttle complete\n')
                fprintf('Transfering aquired data to computer....\n')
                for n = 1:loops
                    fprintf('Start Read %d .... ', n);
                    firstIndex = ((n-1)*numberOfPulses)+1;
                    tic
                    rc = inst.ReadMultipleAdcFrames(chanIndex, firstIndex, numberOfPulses, netArray);
                    assert(rc == 0);
                    samples = double(netArray);
                    fprintf('Read %d Complete\n', n);
                    toc
                    
                    tic
                    %delete mem
                    fprintf('Clear mem %d .... ', n);
                    rc =inst.WipeMultipleAdcFrames(chanIndex, ((n-1)*numberOfPulses)+1, numberOfPulses, 0);
                    assert(rc == 0);
                    fprintf('Clear mem %d Complete\n', n);
                    toc

                    tic 
                    fprintf('Starting interation %d data processing....', n);
                    pulses = reshape(samples, [], numberOfPulses);
                    clear samples;
                    for i = 1:numberOfPulses
                        pulse = pulses(:, i);
%                         pulse = pulse(1024:readLen);
%                         readLen=length(pulse);
                        if i == 209 && n == 1
                            figure(2);clf;
                            plot(pulse);
                            hold on;
                            yline(2048);
                        end  
                        pulseMean = mean(pulse);
                        pulseDC = pulse - pulseMean; % remove DC
                        X = fftshift(fft(pulseDC,padded_len));
                        linFFT = (abs(X)/readLen);
%                         [amp, loc] = max(linFFT);
                        amp=(linFFT(b1) + linFFT(b2));
            %             linFFT_vec(:,i)=linFFT;
                        pulseAmp(i+(numberOfPulses*(n-1))) = amp;
                    end
                    clear pulses;
                    fprintf('Data processing interation %d complete!\n', n);
                    toc
                end

                ivec=1:numberOfPulses*loops;
                time_cycle=pw+(4+2+2+tacq)*1e-6;
                time_axis=time_cycle.*ivec;
                %drop first point
                time_axis(1)=[];pulseAmp(1)=[];
                try
                figure(1);clf;
%                plot(time_axis,(pulseAmp),'b-');hold on;
                plot(time_axis,(smooth(pulseAmp,10)),'r-');
                set(gca,'ylim',[0 max(pulseAmp)*1.1]);  
  
                figure(3);clf;
                plot(pulse);
                hold on;
                yline(2048);
                catch
                    disp('Plot error occured');
                end
                
                %fn=dataBytes; %filename
                a = datestr(now,'yyyy-mm-dd-HHMMSS');
                fn = sprintf([a,'_Proteus']);
                % Save data
                fprintf('Writing data to Z:.....\n');
                save(['Z:\' fn],'pulseAmp','time_axis');
            
            case 4 % Cleanup and save

                rc = inst.SetAdcCaptureEnable(off);
                assert(rc == 0);
                % Free the memory space that was allocated for ADC capture
                rc = inst.FreeAdcReservedSpace();
                assert(rc == 0);
                
%                 clear pulseAmp time_axis;
                fprintf('Complete ready for next aquisition!\n');
            
             case 5 % Shutdown disconnect instr
                 
                 connect = 0;
                 
        end % end switch
        
    end % end while
    
    res = inst.SendScpi(':SYST:ERR?');
    fprintf(1, '\nEnd - server stopped!! \nInstrument Error Status = %s\n', netStrToStr(res.RespStr));
    
    % It is recommended to disconnect from instrumet at the end
    rc = admin.CloseInstrument(instId);    
    
    % Close the administrator at the end ..
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
