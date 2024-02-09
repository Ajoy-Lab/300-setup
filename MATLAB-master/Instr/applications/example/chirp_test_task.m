sampleRateDAC = 9e9;

awg_center_freq = 3.775e9;
                awg_bw_freq = 24e6;
                awg_amp = 1.2;
                sweep_freq = 750;
                sweep_sigma = 0.1;
                symm = 0;
                srs_freq = 0.3625e9; % new value for good chirp
                srs_amp = 3;
                
                  bits = 8;
                
                rampTime = 1/sweep_freq;
                fCenter = awg_center_freq - srs_freq;
                fStart = fCenter - 0.5*awg_bw_freq;
                disp(['fstart = ' num2str(fStart)]);
                fStop = fCenter + 0.5*awg_bw_freq;
                dt = 1/sampleRateDAC;
                dacSignal = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits);   
                fprintf('waveform length - ');
                fprintf(num2str(length(dacSignal)));
                fprintf('\n') ;
                
                t = 0:1/sampleRateDAC:rampTime;
                dacSignal1 = chirp(t,fStart,rampTime,fStop);
                dacSignal2 = fliplr(dacSignal1);
                figure
                pspectrum(dacSignal1,'spectrogram')
                ylim([0.7 0.8])
                figure
                pspectrum(dacSignal2,'spectrogram')
                ylim([0.7 0.8])
                
                
                
function dacWav = makeChirp(sampleRateDAC, rampTime, dt, fStart, fStop, bits)            

    t = 0:1/sampleRateDAC:rampTime;
    dacWave = chirp(t,fStart,rampTime,fStop);
    seglenTrunk = (floor(length(dacWave)/ 64))*64;
    dacWave = dacWave(1:seglenTrunk);
    dacWav = ampScale(bits, dacWave);

end