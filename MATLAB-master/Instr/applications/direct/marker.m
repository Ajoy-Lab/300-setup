% # -*- coding: utf-8 -*-
% """
% Created on Wed Jun 28 13:12:47 2023
% 
% @author: jason
% """

% ### Add trigger marker for function generator ### 

% ## final segment which will contain the marker
chNum = 1;
segNum = 4;

% ## how long the "on" portion needs to be
onLength = 1e-3;

% ## when does the marker start after the last pulse
buffer = 1e-3;

% ## select the final segment
cmd = sprintf(':INST:CHAN %d',chNum);
inst.SendScpi(cmd);
cmd = sprintf(':TRAC:SEL %d',segNum);
inst.SendScpi(cmd);

% ## get the length of the segment
mkrsegment_length = uint8(inst.SendScpi(':TRAC:DEF?'));
mkrsegment_length = uint8(mkrsegment_length/4);

% ## make a new segment
mkr_vector_2 = zeros(mkrsegment_length,1);
mkr_vector_1 = zeros(mkrsegment_length,1);
onLength_points = uint8(onLength*sampleRateDAC/32)*8;
buffer_start = uint8(buffer*sampleRateDAC/32)*8;

% ## make the marker
mkr_vector_on = ones(onLength_points,1);
mkr_vector_2(buffer_start : buffer_start+onLength_points) = mkr_vector_on;
% proteus.inst.timeout = 30000

% # Send the binary-data with *OPC? added to the beginning of its prefix.
mkr_vector = mkr_vector_1 + 2*mkr_vector_2;
% inst.WriteBinaryData('*OPC?; :MARK:DATA', mkr_vector);
inst.WriteBinaryData(':MARK:DATA 0,', mkr_vector);

% % # Set normal timeout
% proteus.inst.timeout = 10000
cmd = ':MARK:SEL 1';
inst.SendScpi(cmd);
cmd = ':MARK:STAT ON';
inst.SendScpi(cmd);
% proteus.checkForError()
