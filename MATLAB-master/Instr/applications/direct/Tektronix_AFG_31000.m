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
            fprintf(obj.gpib_obj, "*CLS");
            fprintf(obj.gpib_obj, "*RST");
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
            fprintf(obj.gpib_obj, sprintf("SOURce1:BURSt:NCYCles %d", ncycles));
            
            % turn on output then set the frequency
            fprintf(obj.gpib_obj, "OUTP1:STAT ON");
            fprintf(obj.gpib_obj, "FUNCTION SIN");
            fprintf(obj.gpib_obj, sprintf("FREQUENCY %d", freq));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:AMPLITUDE %d", Vpp));
            fprintf(obj.gpib_obj, sprintf("VOLTAGE:OFFSET %d", DC_offset));
            fprintf(obj.gpib_obj, sprintf("PHASE:ADJUST %dDEG", phase));
            
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

        function output_off(obj)
            %{
            Turns off the output and resets.
            %}
            fprintf(obj.gpib_obj, "OUTP1:STAT OFF");
            fprintf(obj.gpib_obj, "*CLS");
            fprintf(obj.gpib_obj, "*RST");
        end
        
    end
end