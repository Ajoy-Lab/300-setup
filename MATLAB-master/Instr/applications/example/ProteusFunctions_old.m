classdef ProteusFunctions
    %PROTEUSFUNCTIONS Summary of this class goes here
    %   Detailed explanation goes here

    % Some functions added by Joan Mercade on April 27th 2022
    
    methods

        function [options] = getOptions(~, inst)
            res = inst.SendScpi('*OPT?');  
            options = convertCharsToStrings(char(res.RespStr));
            % options = obj.netStrToStr(res.RespStr);
    
            options = split(options, ',');    
        end

        function [actualChannel] = getActualChannel(~, model_name, labelChannel, cType)
            if contains(model_name, "D")
                if cType == "LAN"
                    actualChannel = labelChannel;
                else
                    labelChannel = mod(labelChannel, 4);
                    switch labelChannel
                        case 1
                            actualChannel = 3;
                        case 2
                            actualChannel = 4;
                        case 3
                            actualChannel = 1;
                        case 4
                            actualChannel = 2;
                    end

                end

            else
                actualChannel = labelChannel;
            end

        end

        function [granul] = getGranularity(~, model_name, options, x1Flag)
            lowGranFlag = false;

            for i = 1:length(options)
                if strncmp(options(i), "G", 1)
                    lowGranFlag = true;
                end
            end
            % Setup model dependant parameters 
            if strncmp(model_name,'P908',4)
                if lowGranFlag
                    granul = 32;
                else
                    granul = 64;
                end
            elseif strncmp(model_name,'P948',4)
                if x1Flag
                    if lowGranFlag
                        granul = 32;
                    else
                        granul = 64;
                    end
                else
                    if lowGranFlag
                        granul = 16;
                    else
                        granul = 32;
                    end
                end
            else
                if lowGranFlag
                    granul = 16;
                else
                    granul = 32;
                end
            end 
        end

        function retval = myQuantization (~,myArray, dacRes, minLevel)  
 
          maxLevel = 2 ^ dacRes - 1;  
          numOfLevels = maxLevel - minLevel + 1;
          
          retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
          retval = retval + minLevel;
          
          retval(retval > maxLevel) = maxLevel;
          retval(retval < minLevel) = minLevel;
        
        end

        function finalWfm = trimGran(~,inWfm, granularity)
            % trimGran - Adjust wfm length for granularity
            %
            % Synopsis
            %   finalWfm = trimGran(inWfm, granularity)
            %
            % Description
            %   Repeat waveform the minmum number of times to meet the
            %   waveform length granularity criteria
            %
            % Inputs ([]s are optional)
            %   (double) inWfm  Input waveform
            %   (int16)  granularity
            %
            % Outputs ([]s are optional)
            %   (double) finalWfm Adjusted waveform
        
            baseL = length(inWfm);
            finaL = lcm(baseL, granularity);
            
            finalWfm = zeros(1, finaL);
            pointer = 1;
            
            while pointer < finaL
                finalWfm(pointer : (pointer + baseL -1)) = inWfm;
                pointer = pointer + baseL;        
            end
        
        end
        
        %% 
        function retval = convertToBynaryOffset (~,myArray, dacRes)
  
          minLevel = 0;
          maxLevel = 2 ^ dacRes - 1;  
          numOfLevels = maxLevel - minLevel + 1;

          retval = round((numOfLevels .* (myArray + 1) - 1) ./ 2);
          retval = retval + minLevel;

          retval(retval > maxLevel) = maxLevel;
          retval(retval < minLevel) = minLevel;

        end
        %%
        function [normI,  normQ] = NormalIq(~,wfmI, wfmQ)    
    
            maxPwr = max(wfmI.*wfmI + wfmQ .* wfmQ);
            maxPwr = maxPwr ^ 0.5;

            normI = wfmI ./ maxPwr;
            normQ = wfmQ ./ maxPwr;

        end
        %% 
        function [fact] = GetNormalIqFactor(~,wfmI, wfmQ)    
    
            maxPwr = max(wfmI.*wfmI + wfmQ .* wfmQ);
            fact = maxPwr ^ 0.5;    

        end
        %% 
        function outWfm = Interleave(~,wfmI, wfmQ)   

            wfmLength = length(wfmI);
            outWfm = zeros(1, 2 * wfmLength);

            outWfm(1:2:(2 * wfmLength - 1)) = wfmI;
            outWfm(2:2:(2 * wfmLength)) = wfmQ;
        end
        %% 
        function dataOut = expanData(~,inputWfm, oversampling)
            dataOut = zeros(1, oversampling * length(inputWfm));
            dataOut(1:oversampling:length(dataOut)) = inputWfm;
        end
        
        %% 
        function dataOut = getRnData(~,nOfS, bPerS)
            maxVal = 2 ^ bPerS;
            dataOut = maxVal * rand(1, nOfS);
            dataOut = floor(dataOut);
            dataOut(dataOut >= maxVal) = maxVal - 1;    
        end
        %% 

        function retval = MyIdealInterpolation2 (~,myArray, xFactor, quality, bwFraction)
          % Function used in traditional resampling
          
          %expansion by zero-padding
          retval = zeros(1, xFactor * length(myArray));
          retval([1:xFactor:end]) = myArray;
          % "Ideal" Interpolation filter
          lenSinc = quality; 
          
          mySinc = -lenSinc * xFactor : 1 : lenSinc * xFactor ;  
          mySinc = sinc(bwFraction .* mySinc / xFactor);
          %The FIR taps must be normalized so the DC response is 0dB
          normFactor = sum(mySinc(1:xFactor:length(mySinc)));
          normFactor = 1/ normFactor;
          mySinc = normFactor .* mySinc;
          myWindow = blackman(length(mySinc));
          myWindow = myWindow.'; 
          mySinc = mySinc .* myWindow;
          %convolution
          retval = cconv(retval, mySinc, length(retval));
          %retval = real(retval);
          
        end

        function retval = MyOldResampleWfm(obj, myArray, quality, inSr, outSr)

            % Traditional resampling function. Just for reference and comparison
        
            [inSr, outSr] = obj.reduceFraction(inSr, outSr);
            
            if outSr > inSr
                bwFraction = 1.0;
            else
                bwFraction = outSr / inSr;
            end
            
            bwFraction = 0.98 * bwFraction;  
            retval = obj.MyIdealInterpolation2 (myArray, outSr, quality, bwFraction);    
            retval = retval(1:inSr:length(retval));   
        end

        function waveform = File_Modulation_wave(obj,modType,fileI, fileQ,fileFormat,sampleRateIn,sampleRateOut,interpol)
           % Read data from file
           % modType -2 = One File IQ
           % modType -3 = Two Files
           % fileFormat "mat" "csv"

           quality = 100;
           waveform = double.empty;
           sampleRateOut = sampleRateOut / interpol;

           [inSr, outSr] = obj.reduceFraction(sampleRateIn, sampleRateOut);

           if modType == -3
               if fileFormat == "csv"
                   wfmI = readmatrix(fileI);
                   wfmI = wfmI';
                   wfmQ = readmatrix(fileQ);
                   wfmQ = wfmQ';
                   
               elseif fileFormat == "mat"

               end
           elseif modType == -2
               if fileFormat == "csv"
                   wfmI = readmatrix(fileI);
                   wfmI = wfmI';
                   wfmQ = hilbert(wfmI);
                   wfmI = real(wfmQ);
                   wfmQ = imag(wfmQ);  

               elseif fileFormat == "mat"
               end
           end

           wfmI = obj.MyOldResampleWfm(wfmI, quality, sampleRateIn, sampleRateOut);
           wfmQ = obj.MyOldResampleWfm(wfmQ, quality, sampleRateIn, sampleRateOut);                   
           waveform = wfmI + 1i * wfmQ;

        end

        function waveform = Analog_Modulation_wave(obj,modType,minCycles,sampleRate,interpol,granul,param1,param2)
   
            %ANALOG MODULATION WAVEFORM CALCULATION
            % modType = 0, AM; 1, FM; 2, PM; 3, SSB; 4, CHIRP;    

            %AM SETTINGS
            amModIndex = param1; %Modulation Index in %
            amModFreq = param2; %Modulation frequency in HZ

            %FM SETTINGS
            fmFreqDev = param1; %Peak Frequency Deviaton in Hz
            fmModFreq = param2; %Modulaition frequency in HZ

            %PM SETTINGS
            pmPhaseDev = param1; %Peak Phase Deviaton in Rads
            pmModFreq = param2; %Modulation frequency in HZ

            %SSB SETTINGS
            ssbModFreq = param2; %Modulation frequency in HZ

            %CHIRP SETTINGS
            chirpSweepRange = param1;
            chirpSweepTime = param2;

            %Waveform Length Calculation    
            modFreq = amModFreq;

            if modType == 1
                modFreq = fmModFreq;
            elseif modType == 2
                modFreq = pmModFreq;
            elseif modType == 3
                modFreq = ssbModFreq;
            elseif modType == 4
                modFreq = 1 / chirpSweepTime;
            end

            actualSR = sampleRate / interpol;
            if modType ~= 4 
                numOfSamples = round(actualSR / abs(modFreq / minCycles));
            else
                numOfSamples = round(actualSR / abs(modFreq));
            end
            totalNumOfSamples = numOfSamples;

            % As samples sent to the instrument are twice the number of complex
            % samples, granul must be defined as half the actual number

            while modType ~= 4 && mod(totalNumOfSamples, granul) ~= 0
                totalNumOfSamples = totalNumOfSamples + numOfSamples;
            end

            numOfSamples = totalNumOfSamples;
            fRes = actualSR / numOfSamples;

            % Round modFreq to the nearest integer number of Cycles

            modFreq = round(modFreq / fRes) * fRes;

            %Waveform calculation
            fprintf(1, 'WAVEFORM CALCULATION\n');

            waveform = 0: (numOfSamples - 1);
            waveform = (1 / actualSR) .* waveform;
            waveform = waveform - (numOfSamples / (2 * actualSR));

            if modType == 0
                waveform = 1 + amModIndex/100 .* sin(2 * pi * modFreq * waveform);
                waveform = waveform + 1i * waveform;
            elseif modType == 1
                fmFreqDev = round(fmFreqDev / fRes) * fRes;
                freqInst = fmFreqDev / modFreq * sin(2 * pi * modFreq * waveform);
                waveform = cos(freqInst) + 1i * sin(freqInst); 
                clear freqInst;
            elseif modType == 2
                phaseInst = pmPhaseDev * sin(2 * pi * modFreq * waveform);
                waveform = cos(phaseInst) + 1i * sin(phaseInst);
                clear phaseInst;
            elseif modType == 3
                waveform = 2 * pi * modFreq * waveform;
                waveform = cos(waveform) + 1i * sin(waveform);
            elseif modType == 4
                chirpSweepRange = chirpSweepRange / 2;
                chirpSweepRange = round(chirpSweepRange / fRes) * fRes;
                freqInst = (actualSR * chirpSweepRange / numOfSamples) * waveform;
                freqInst = 2 * pi * freqInst .* waveform;
                waveform = cos(freqInst) + 1i * sin(freqInst);
                clear freqInst;    
                waveform = obj.trimGran(waveform, granul);    
            end

            % waveform conditioning:    
            waveform = waveform./((mean(abs(waveform).^2))^0.5);

        end

        function retval = IqModulator(  ~, ...
                                        iqWfm, ...
                                        cFreq, ...
                                        cPhase, ...
                                        sRate)
  
          % The carrier complex waveform is calculated
          carrier = 0:(length(iqWfm) - 1);
          carrier = carrier / sRate;
          %
          %carrier = cos(2*pi*cFreq * carrier) + 1i * sin(2*pi*cFreq * carrier);
          carrier = exp(1i * 2 * pi * cFreq * carrier +  cPhase);
          % IQ modulation is performed by complex multiplication of the complex baseband
          % and the complex carrier
          retval = iqWfm .* carrier;
          % The actual modulated waveform is just the real part of the product
          retval = real(retval);
        
        end
        
        %% 
        function [dataOut] = Digital_Modulation_wave(obj,modType,numOfSymbols,symbolRate,rollOff,sampleRate,interpol)      
            % modType     Modulation
            % 5                 QPSK
            % 6                 QAM16
            % 7                 QAM32
            % 8                 QAM64
            % 9                 QAM128
            %10                 QAM256
            %11                 QAM512
            %12                 QAM1024
            
            if modType == 5
                bitsPerSymbol = 2;
            elseif modType == 6
                bitsPerSymbol = 4;
            elseif modType == 7
                bitsPerSymbol = 5;
            elseif modType == 8
                bitsPerSymbol = 6;
            elseif modType == 9
                bitsPerSymbol = 7;
            elseif modType == 10
                bitsPerSymbol = 8;
            elseif modType == 11
                bitsPerSymbol = 9;
            elseif modType == 12
                bitsPerSymbol = 10;
            else
                bitsPerSymbol = 2;
            end
       

            % Waveform Length Calculation
            sampleRate = sampleRate / interpol;

            [decimation, oversampling] = obj.reduceFraction(symbolRate, sampleRate);              


            % Create IQ for QPSK/QAM    
            % accuracy is the length of-1 the shaping filter
            accuracy = 64;
            fType = 'normal'; % 'normal' or 'sqrt'
            % Get symbols in the range 1..2^bps
            data = obj.getRnData(numOfSymbols, bitsPerSymbol);
            % Map symbols to I/Q constellation locations
            [dataI, dataQ] = obj.getIqMap(data, bitsPerSymbol);
            % Adapt I/Q sample rate to the AWG's

            dataI = obj.expanData(dataI, oversampling);
            dataQ = obj.expanData(dataQ, oversampling);
            % Calculate baseband shaping filter
            rsFilter = rcosdesign(rollOff,accuracy,oversampling, fType);
            % Apply filter through circular convolution
            dataI = cconv(dataI, rsFilter, length(dataI));
            dataQ = cconv(dataQ, rsFilter, length(dataQ));

            dataI = dataI(1:decimation:length(dataI));
            dataQ = dataQ(1:decimation:length(dataQ));
            
            L = length(dataI);
            Lr = 32*floor(L/32);
            
            dataI = dataI(1:Lr);
            dataQ = dataQ(1:Lr);
            
            
            % Output waveforfm must be made of complex samples
            dataOut = dataI + 1i * dataQ;
        end
        %% 
        function [outNum, outDen] = reduceFraction(~,num, den)
        %reduceFraction Reduce num/den fraction
        %   Use integers although not mandatory
            num = round(num);
            den = round(den);
            % Reduction is obtained by calcultaing the greater common divider...
            G = gcd(num, den);
            % ... and then dividing num and den by it.
            outNum = num / G;
            outDen = den / G;
        end
        %% 
        function [symbI, symbQ] = getIqMap(~,data, bPerS)
   
            if bPerS == 5 % QAM32 mapping
                lev = 6;
                data = data + 1;
                data(data > 4) = data(data > 4) + 1;
                data(data > 29) = data(data > 29) + 1; 

            elseif bPerS == 7 % QAM128 mapping      
                lev = 12;
                data = data + 2;
                data(data > 9) = data(data > 9) + 4;
                data(data > 21) = data(data > 21) + 2;
                data(data > 119) = data(data > 119) + 2;
                data(data > 129) = data(data > 129) + 4;

             elseif bPerS == 9 % QAM512 mapping       
                lev = 24;
                data = data + 4;
                data(data > 19) = data(data > 19) + 8;
                data(data > 43) = data(data > 43) + 8;
                data(data > 67) = data(data > 67) + 8;
                data(data > 91) = data(data > 91) + 4;
                data(data > 479) = data(data > 479) + 4;
                data(data > 499) = data(data > 499) + 8;
                data(data > 523) = data(data > 523) + 8;
                data(data > 547) = data(data > 547) + 8;            
            else
                lev = 2 ^ (bPerS / 2); % QPSK, QAM16, QAM64, QAM256, QAM1024      
            end

            symbI = floor(data / lev);
            symbQ = mod(data, lev);
            lev = lev / 2 - 0.5;   
            symbI = (symbI - lev) / lev;
            symbQ = (symbQ - lev) / lev;
        end
        
        %%
        function [str] = netStrToStr(~,netStr)
            try
                str = convertCharsToStrings(char(netStr));
            catch        
                str = '';
            end
        end

        %% 
        function [dataout] = ConvertSampleToNormalsigned(~,inp,size)
            M = 2^(size-1);
            
            dataout = double(inp) - M;
            dataout = dataout ./ M;
        end
        
        %% 
        function [env,gauss_i,gauss_q] = gauss_env(~,pw,pl,fs,fc,interp,phase,direct,direct_lo,SQP,NP,PG)
    
            res = 64;
           
            fs = fs / interp;
            sigma = double(pw / 6);
            variance = sigma^2;
            pg = PG;
            wavelength = pl * fs;
            wavelength = res * ceil(wavelength / res);
            N = wavelength;
            ts = 1 / fs;
            t = linspace(-N*ts/2, (N-2)*ts/2, N);
            
            
            phase = phase * pi / 180;
            if NP > 1
                fc_v = linspace(fc, fc+NP*pg, NP);
                fc_v = fc_v(:);     % convert to coloumn
            else
                fc_v = fc;
            end
            
            M = t.*fc_v;
            sinWave_m = sin(phase + 2 * pi * M);
            cosWave_m = cos(phase + 2 * pi * M); 

            gauss_sq_pulse = zeros(1,N);
            gauss_sq_pulse(1:int32(pw/ts)) = 1;
            gauss_e = exp(-t.^2/(2*variance));

            
            if SQP == false
                gauss_i_m = cosWave_m .* gauss_e;
                gauss_q_m = sinWave_m .* gauss_e;
            else
                gauss_i_m = cosWave_m .* gauss_sq_pulse;
                gauss_q_m = sinWave_m .* gauss_sq_pulse;
            end

            flo = direct_lo;
            lo_sinWave = sin(2 * pi * flo .* t);
            lo_cosWave = cos(2 * pi * flo .* t);

            mod_m = gauss_i_m .* lo_cosWave - gauss_q_m .* lo_sinWave;
            
            
            if NP > 1
                mod = sum(mod_m);
            else
                mod = mod_m;
            end

            if direct == true
                env = mod; 
            else
                if SQP == false
                    env = gauss_e;
                else
                    env = gauss_sq_pulse;
                end
            end

            if NP > 1
                gauss_i = sum(gauss_i_m);
                gauss_q = sum(gauss_q_m);
            else
                gauss_i = gauss_i_m;
                gauss_q = gauss_q_m;
            end
           
        end
        
        %% 
        function [out] = clip(~,in,min,max)
            in(in>max) = max;
            in(in<min) = min;
            out = in;
        end
        
        %% 
        function [k_i,k_q] =  iq_kernel(~,fs,flo,coe_file_path)
            % load coe data for the FIR filter
            data = readtable(coe_file_path);
            coe = table2array(data);
            TAP = length(coe);
            fprintf('loaded %d TAP filter from %s',TAP,coe_file_path);
            L = 10240;
            k = ones(1,L+TAP);

            ts = 1 / fs;
            t = linspace(0, (L-1)*ts, L);

            loi = cos(2 .* pi .* flo .* t);
            loq = -(sin(2 .* pi .* flo .* t));

            k_i = zeros(1, L);
            k_q = zeros(1, L);

            for l = 1:L
                b = 0;
                for n = 1:TAP
                    b = b + k(l+n)*coe(n);
                end
                k_q(l) = loq(l) * b;
                k_i(l) = loi(l) * b;
            end

        end
        
        %% 
        function [out_i , out_q] =  convertIQtoFixSignedSample(obj,inp_i,inp_q,size)
            out_i = zeros(1,length(inp_i));
            out_i = typecast(out_i,'uint32');
            
            out_q = zeros(1,length(inp_q));
            out_q = typecast(out_q,'uint32');

            [inp_i,inp_q] = obj.NormalIq(inp_i,inp_q);

            M = 2^(size-1);
            A = 2^(size);

            for i = 1:length(inp_i)
                if inp_i(i) < 0
                    out_i(i) = round(inp_i(i)*M) + A;
                else
                    out_i(i) = round(inp_i(i)*(M-1));
                end
            end

            for i = 1:length(inp_q)
                if inp_q(i) < 0
                    out_q(i) = round(inp_q(i)*M) + A;
                else
                    out_q(i) = round(inp_q(i)*(M-1));
                end
            end            
        end        
        
        %% 
        function hex_str = Convert_FIX_m_n_to_hex(~,val,fp_size)
            num = val * (2^fp_size);
            num = int32(num);
            hex_str = dec2hex(num);
            hex_str = strcat('0x',hex_str);    
        end
        
        %% 
        function [kernel_data] = pack_kernel_data(obj,ki,kq,EXPORT)
            L = length(ki)/5;
            out_i = zeros(1,L*4);
            out_q = zeros(1,L*4);
            
            % convert the signed number into 12bit FIX1_11 presentation
            [b_ki,b_kq] = obj.convertIQtoFixSignedSample(ki,kq,12);            

            % convert 12bit to 15bit because of FPGA memory structure
            for i = 1:L
                j = i - 1;
                s1 = bitand(b_ki(j*5+2),7	) * 4096 + bitand(b_ki(j*5+1),4095)         ;             
                s2 = bitand(b_ki(j*5+3),63	) *  512 + bitand(b_ki(j*5+2),4088) / 2^3   ; 
                s3 = bitand(b_ki(j*5+4),511	) *   64 + bitand(b_ki(j*5+3),4032) / 2^6   ;
                s4 = bitand(b_ki(j*5+5),4095) *    8 + bitand(b_ki(j*5+4),3584) / 2^9   ; 
                out_i(1+j*4) = s1;
                out_i(2+j*4) = s2;
                out_i(3+j*4) = s3;
                out_i(4+j*4) = s4;
            end
            
            for i = 1:L
                j = i - 1;
                s1 = bitand(b_kq(j*5+2),7	) * 4096 + bitand(b_kq(j*5+1),4095)         ;               
                s2 = bitand(b_kq(j*5+3),63	) *  512 + bitand(b_kq(j*5+2),4088) / 2^3   ;
                s3 = bitand(b_kq(j*5+4),511	) *   64 + bitand(b_kq(j*5+3),4032) / 2^6   ;
                s4 = bitand(b_kq(j*5+5),4095) *    8 + bitand(b_kq(j*5+4),3584) / 2^9   ;
                out_q(1+j*4) = s1;
                out_q(2+j*4) = s2;
                out_q(3+j*4) = s3;
                out_q(4+j*4) = s4;
            end            

            fout_i = zeros(1,L*4);
            fout_q = zeros(1,L*4);

            for i = 1:(L*4)
                if out_i(i) >16383
                    fout_i(i) = out_i(i) - 32768;
                else
                    fout_i(i) = out_i(i);
                end
            end
            
            for i = 1:(L*4)
                if out_q(i) >16383
                    fout_q(i) = out_q(i) - 32768;
                else
                    fout_q(i) = out_q(i);
                end
            end            

            kernel_data = out_q.*(2^16) + out_i;
           
            
            kernel_filt_data = [fout_i;fout_q];
            kernel_filt_data = kernel_filt_data';

            sim_kernel_data = compose("%X",kernel_data);
            
            if EXPORT == true
                [~, ~, ~] = mkdir('kernel_output_files');
                writematrix(kernel_filt_data,".\kernel_output_files\kernel_filt.csv");
                writematrix(kernel_data',".\kernel_output_files\mem_data.csv");
                writematrix(sim_kernel_data',".\kernel_output_files\sim_mem_data.csv");
            end            
        end                  
    end
end


    
  