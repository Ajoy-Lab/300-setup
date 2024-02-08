function [myWaveI, myWaveQ] = makeSqPulse(modFreq, pulseLen, amplitude, phase, mods, sampleRateDAC)
        
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