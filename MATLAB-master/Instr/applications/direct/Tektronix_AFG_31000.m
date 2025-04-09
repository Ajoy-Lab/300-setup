classdef Tektronix_AFG_31000
    properties
        gpib_obj
        port_name
    end
    methods
        function obj=Tektronix_AFG_31000(port_name)
            %{ 
            initialize Tektrnoix (AFG 31000 Arbitrary function generator) 
            device to generate arbitrary AC field for AC sensing 
            %}
            obj.gpib_obj=visa('ni',port_name);
            fopen(obj.gpib_obj);
            %fprintf(obj.gpib_obj, "*CLS");
            %fprintf(obj.gpib_obj, "*RST");
        end
        
        function sweep(obj, fstart, fstop, sweep_t, Vpp)
             % set souce1
            fprintf(obj.gpib_obj, "SOURce1:FREQuency:MODE SWE");
            fprintf(obj.gpib_obj, sprintf("SOUR1:FREQ:STAR %.3f", fstart));
            fprintf(obj.gpib_obj, sprintf("SOUR1:FREQ:STOP %.3f", fstop));
            fprintf(obj.gpib_obj, sprintf("SOUR1:SWE:TIME %.3f", sweep_t));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT:LEV:IMM:AMP %.3f", Vpp));
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
            fprintf(obj.gpib_obj, "TRIG:SLOP POS");
            fprintf(obj.gpib_obj, "SOUR1:SWEep:MODE MAN");
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
        end
        
        function init_AFG_RF(obj, freq, Vpp, DC_offset, phase)
            fprintf(obj.gpib_obj, "SOUR1:FUNC SIN");
            fprintf(obj.gpib_obj, sprintf("SOUR1:FREQ %.3f", freq));
            %disp(sprintf("SOUR1:FREQ %.3f", freq));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT %.3f", Vpp));
            %disp(sprintf("SOUR1:VOLT %.3f", Vpp));
            fprintf(obj.gpib_obj, sprintf("SOUR1:PHAS %dDEG", phase));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT:LEV:IMM:OFFS %dV", DC_offset));
        end
        
        function init_AFG_RF_FM(obj, freq, Vpp, DC_offset, phase, FMshape, FMfreq, FMdeviation, external_modulation)
            fprintf(obj.gpib_obj, "SOUR1:FUNC SIN");
            fprintf(obj.gpib_obj, "SOUR1:FM:STAT ON");
            
            fprintf(obj.gpib_obj, sprintf("SOUR1:FREQ %.3f", freq));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT %.3f", Vpp));
            fprintf(obj.gpib_obj, sprintf("SOUR1:PHAS %dDEG", phase));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT:LEV:IMM:OFFS %dV", DC_offset));
            
            fprintf(obj.gpib_obj, "SOUR1:FM:SOUR INT");
            fprintf(obj.gpib_obj, sprintf("SOUR1:FM:INT:FUNC %s", FMshape));
            fprintf(obj.gpib_obj, sprintf("SOUR1:FM:INT:FREQ %.3f", FMfreq));
            disp(sprintf("SOUR1:FM:DEV %.3fHz", FMdeviation));
            fprintf(obj.gpib_obj, sprintf("SOUR1:FM:DEV %.3fHz", FMdeviation));
            
            
            if external_modulation
                fprintf(obj.gpib_obj, "SOUR1:FM:SOUR EXT");
                %fprintf(obj.gpib_obj, "SOUR2:FUNC SIN");    % noise function PRN
                %fprintf(obj.gpib_obj, "SOUR2:VOLT:LEV:IMM:AMPL 1Vpp");
            end
            
        end
        
        function init_AFG_RF_twoFM(obj, freq, Vpp, DC_offset, phase, FMshape, FMfreq, FMdeviation, external_modulation)
            fprintf(obj.gpib_obj, "SOUR1:FUNC SIN");
            fprintf(obj.gpib_obj, "SOUR2:FUNC SIN");
            fprintf(obj.gpib_obj, "SOUR1:FM:STAT ON");
            fprintf(obj.gpib_obj, "SOUR2:FM:STAT ON");
            
            freq_ch1 = freq+FMdeviation;
            freq_ch2 = freq+3*FMdeviation;
            fprintf(obj.gpib_obj, sprintf("SOUR1:FREQ %.3f", freq_ch1));
            fprintf(obj.gpib_obj, sprintf("SOUR2:FREQ %.3f", freq_ch2));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT %.3f", Vpp));
            fprintf(obj.gpib_obj, sprintf("SOUR2:VOLT %.3f", Vpp));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT:OFFS %dV", DC_offset));
            fprintf(obj.gpib_obj, sprintf("SOUR2:VOLT:OFFS %dV", DC_offset));
            
            fprintf(obj.gpib_obj, "SOUR1:FM:SOUR INT");
            fprintf(obj.gpib_obj, "SOUR2:FM:SOUR INT");
            fprintf(obj.gpib_obj, sprintf("SOUR1:FM:INT:FUNC %s", FMshape));
            fprintf(obj.gpib_obj, sprintf("SOUR2:FM:INT:FUNC %s", FMshape));
            fprintf(obj.gpib_obj, sprintf("SOUR1:FM:INT:FREQ %.3f", FMfreq));
            fprintf(obj.gpib_obj, sprintf("SOUR2:FM:INT:FREQ %.3f", FMfreq*1.07));
            disp(sprintf("SOUR1:FM:DEV %.3fHz", FMdeviation));
            fprintf(obj.gpib_obj, sprintf("SOUR1:FM:DEV %.3fHz", FMdeviation));
            fprintf(obj.gpib_obj, sprintf("SOUR2:FM:DEV %.3fHz", FMdeviation));
            
            if external_modulation
                fprintf(obj.gpib_obj, "SOUR1:FM:SOUR EXT");
                %fprintf(obj.gpib_obj, "SOUR2:FUNC SIN");    % noise function PRN
                %fprintf(obj.gpib_obj, "SOUR2:VOLT:LEV:IMM:AMPL 1Vpp");
            end
            
        end


        function init_AFG_RF_AM(obj, freq, Vpp, DC_offset, phase)
            fprintf(obj.gpib_obj, "SOUR1:FUNC SIN");
%             fprintf(obj.gpib_obj, "SOUR1:FM:STAT ON");
            
            fprintf(obj.gpib_obj, sprintf("SOUR1:FREQ %.3f", freq));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT %.3f", Vpp));
            fprintf(obj.gpib_obj, sprintf("SOUR1:PHAS %dDEG", phase));
            fprintf(obj.gpib_obj, sprintf("SOUR1:VOLT:LEV:IMM:OFFS %dV", DC_offset));
            
            fprintf(obj.gpib_obj, "SOUR1:AM:STAT ON");
            fprintf(obj.gpib_obj, "SOUR1:AM:SOUR EXT");
            fprintf(obj.gpib_obj, "SOUR1:AM:DEPT 120");
            
        
            
        end
        
        function burst_mode_trig_sinwave(obj, freq, Vpp, DC_offset, phase, ncycles, add_external)
            %{
            1. Sets Tektronix device to burst mode, waiting for a trigger
            2. Once trigger comes, output a sine wave with given Vpp, freq,
               phase, DCoffset and phase
            %}
            % set trigger souce

            if nargin < 7
                add_external = false;
            end
            
            if add_external
                fprintf(obj.gpib_obj, 'SOURce1:COMBine:FEED "EXTernal"');
            end 
            
            fprintf(obj.gpib_obj, "TRIG:SLOP POS");
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
            
            % set souce1
            fprintf(obj.gpib_obj, "BURSt:STATE ON");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:INFInite:REARm");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:IDLE DC");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:MODE TRIG");
            % for now set it as infinite cycles
%             if ncycles == "INF"
%                 fprintf(obj.gpib_obj, "SOURce1:BURSt:NCYCles INF");
%             else
%                 fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", ncycles));
%             end
            fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", ncycles));
            
            % turn on output then set the frequency
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
            fprintf(obj.gpib_obj, "FUNCTION SIN");
            fprintf(obj.gpib_obj, sprintf("FREQUENCY %d", freq));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:AMPLITUDE %d", Vpp));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:OFFSET %d", DC_offset));
            fprintf(obj.gpib_obj, sprintf("PHASE:ADJUST %dDEG", phase));

            
        end
        
        function burst_mode_trig_rectwave(obj, freq, Vpp, DC_offset, phase, ncycles, add_external)
            %{
            1. Sets Tektronix device to burst mode, waiting for a trigger
            2. Once trigger comes, output a sine wave with given Vpp, freq,
               phase, DCoffset and phase
            %}
            % set trigger souce

            if nargin < 7
                add_external = false;
            end
            
            if add_external
                fprintf(obj.gpib_obj, 'SOURce1:COMBine:FEED "EXTernal"');
            end 
            
            fprintf(obj.gpib_obj, "TRIG:SLOP POS");
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
            
            %%% set run mode to triggered
            %%fprintf(obj.gpib_obj, "SEQControl:RMODe TRIGgered");
            % stop output and wait for the next trigger
            % fprintf(obj.gpib_obj, "SOURce1:BURSt:INFInite:REARm");
            
            % set souce1
            fprintf(obj.gpib_obj, "BURSt:STATE ON");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:INFInite:REARm");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:IDLE DC");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:MODE TRIG");
            % for now set it as infinite cycles
%             if ncycles == "INF"
%                 fprintf(obj.gpib_obj, "SOURce1:BURSt:NCYCles INF");
%             else
%                 fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", ncycles));
%             end
            fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", ncycles));
            
            % turn on output then set the frequency
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
            fprintf(obj.gpib_obj, "FUNCTION SQU");
            fprintf(obj.gpib_obj, sprintf("FREQUENCY %d", freq));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:AMPLITUDE %d", Vpp));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:OFFSET %d", DC_offset));
            fprintf(obj.gpib_obj, sprintf("PHASE:ADJUST %dDEG", phase));
                        
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR TIM");
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
            
        end

        function burst_mode_trig_waveform(obj, waveform, freq, Vpp, DC_offset, phase, ncycles, add_external)
            %{
            1. Sets Tektronix device to burst mode, waiting for a trigger
            2. Once trigger comes, output the specified waveform (sine or square) with given Vpp, freq, phase, DC_offset, and phase
            %}
            
            % Set default value for add_external if not provided
            if nargin < 8
                add_external = false;
            end
            
            % Set trigger source
            if add_external
                fprintf(obj.gpib_obj, 'SOURce1:COMBine:FEED "EXTernal"');
            end
            
            fprintf(obj.gpib_obj, "TRIG:SLOP POS");
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
            
            % Set source1
            fprintf(obj.gpib_obj, "BURSt:STATE ON");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:INFInite:REARm");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:IDLE DC");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:MODE TRIG");
            fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", ncycles));
            
            % Turn on output then set the frequency and waveform
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
            
            % Use switch-case to set the waveform
            switch waveform
                case 'SIN'
                    fprintf(obj.gpib_obj, "FUNCTION SIN");
                case 'SQU'
                    fprintf(obj.gpib_obj, "FUNCTION SQU");
                case 'TRI'
                    fprintf(obj.gpib_obj, "FUNCTION RAMP");
                otherwise
                    error('Unsupported waveform type. Use "SIN" for sine wave or "SQU" for square wave.');
            end
            
            fprintf(obj.gpib_obj, sprintf("FREQUENCY %d", freq));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:AMPLITUDE %d", Vpp));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:OFFSET %f", DC_offset));
            fprintf(obj.gpib_obj, sprintf("PHASE:ADJUST %dDEG", phase));
            
            % Additional code to activate trigger (don't know why it is necessary but it is)
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR TIM");
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
        end
        
        function burst_mode_DC(obj, voltage)
            %{
            1. Sets Tektronix device to burst mode, waiting for a trigger
            2. Once trigger comes, output a effective DC field with voltage
            as voltage
            %}
            % set trigger souce
            fprintf(obj.gpib_obj, "TRIG:SLOP POS");
            fprintf(obj.gpib_obj, "TRIG:SEQ:SOUR EXT");
            
            % set souce1
            fprintf(obj.gpib_obj, "BURSt:STATE ON");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:IDLE DC");
            fprintf(obj.gpib_obj, "SOURce1:BURSt:MODE GATE");
            % for now set it as infinite cycles
            fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", 100));
            
            % turn on output then set the frequency
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
            fprintf(obj.gpib_obj, "FUNCTION SQU");
            fprintf(obj.gpib_obj, sprintf("FREQUENCY %d", 100));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:AMPLITUDE %d", voltage));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:OFFSET %d", 0));
            fprintf(obj.gpib_obj, sprintf("PHASE:ADJUST %dDEG", 0));
            fprintf(obj.gpib_obj, "SOURce1:BURSt:INFInite:REARm");
        end
        
        function apply_DC(obj, voltage)
            fprintf(obj.gpib_obj, "SOURce1:FUNCtion:SHAPe DC");
            fprintf(obj.gpib_obj, sprintf("SOURce1:VOLT:LEV:IMM:OFFS %d", voltage));
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
        end
    
        function output_off(obj)
            %{
            Turns off the output and resets.
            %}
            fprintf(obj.gpib_obj, "OUTP1:STAT OFF");
            fprintf(obj.gpib_obj, "*CLS");
            fprintf(obj.gpib_obj, "*RST");
        end
        
        function just_output_off(obj)
            %{
            Turns off the output.
            %}
            fprintf(obj.gpib_obj, "OUTP1:STAT OFF");
        end
        
    end
end