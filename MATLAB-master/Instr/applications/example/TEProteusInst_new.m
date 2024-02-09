% =========================================================================
% Copyright (C) 2016 Tabor-Electronics Ltd <http://www.taborelec.com/>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>. 
% =========================================================================
% Author: Nadav Manos
% Date: Nov 23, 2016
% Version: 1.0.1
% $Revision: 3336 $

classdef TEProteusInst < handle
    % TEProteusInst: NI-VISA based connection to Proteus Instrument.
    
   
    properties
        ParanoiaLevel = 1; % Paranoia level (0:low, 1:normal, 2:high)        
    end
    
    properties (SetAccess=private)
        ConnStr = ''; % The Connection-String        
        ViSessn = 0;  % VISA Session
    end
    
    properties (Constant=true)
        VISA_IN_BUFF_SIZE = 8192;   % VISA Input-Buffer Size (bytes)
        VISA_OUT_BUFF_SIZE = 8192;  % VISA Output-Buffer Size (bytes)
        VISA_TIMEOUT_SECONDS = 10;  % VISA Timeout (seconds)
        BINARY_CHUNK_SIZE = 8192;   % Binary-Data Write Chunk Size (bytes)
        WAIT_PAUSE_SEC = 0.02;      % Waiting pause (seconds)
    end
    
    methods % public
        
        function obj = TEProteusInst(connStr, paranoiaLevel)
            % TEProteusInst - Handle Class Constructor
            %
            % Synopsis
            %   obj = TEProteusInst(connStr, [verifyLevel])
            %
            % Description
            %   This is the constructor of the VisaConn (handle) class.
            %
            % Inputs ([]s are optional)
            %   (string) connStr      connection string: either a full  
            %                         VISA resource name, or an IP-Address.
            %   (int) [paranoiaLevel = 1] paranoia level [0,1 or 2].
            % 
            % Outputs
            %   (class) obj      VisaConn class (handle) object.
            %
            
            assert(nargin == 1 || nargin == 2);
            
            ipv4 = '^(?:[0-9]{1,3}\.){3}[0-9]{1,3}$';
            if 1 == regexp(connStr, ipv4)
                connStr = sprintf('TCPIP0::%s::5025::SOCKET', connStr);
            end
            
            if nargin == 2
                %verifyLevel = varargin(1);
                if paranoiaLevel < 1
                    obj.ParanoiaLevel = 0;
                elseif paranoiaLevel > 2
                    obj.ParanoiaLevel = 2;
                else
                    obj.ParanoiaLevel = fix(paranoiaLevel);
                end
            else
                obj.ParanoiaLevel = 1;
            end
            
            obj.ConnStr = connStr;
            
            obj.ViSessn = visa('NI', connStr);
            set(obj.ViSessn, 'OutputBufferSize', obj.VISA_OUT_BUFF_SIZE);
            set(obj.ViSessn, 'InputBufferSize', obj.VISA_IN_BUFF_SIZE);
            obj.ViSessn.Timeout = obj.VISA_TIMEOUT_SECONDS;
            %obj.ViSessn.Terminator = newline;
           
        end
        
        function delete(obj)
            % delete - Handle Class Destructor
            %
            % Synopsis
            %   obj.delete()
            %
            % Description
            %   This is the destructor of the VisaConn (handle) class.
            %   (to be called on a VisaConn class object).
            %
            % Examples
            %   connStr = 'TCPIP::192.168.0.170::5025::SOCKET';
            %   inst = TEWXAwg(connStr, HighLevel);
            %   inst.Connect();
            %
            %   % do some work with the instrument ...
            %
            %   inst.delete(); % delete it when it is no more needed.
            
            obj.Disconnect();
            delete(obj.ViSessn);
            obj.ViSessn = 0;
        end
        
        function ok = Connect(obj)
            % Connect - open connection to remote instrument.
            %
            % Synopsis
            %    ok = obj.Connect()
            %
            % Description
            %    Open connection to the remote instrument
            %
            % Outputs
            %    (boolean) ok   true if succeeded; otherwise false.
            %
            % Examples
            %   connStr = 'TCPIP::192.168.0.170::5025::SOCKET';
            %   inst = TEWXAwg(connStr, HighLevel);
            %   inst.Connect();
            %
            %   % do some work with the instrument ...
            
            ok = false;
            try
                if strcmp(obj.ViSessn.Status, 'open')
                    ok = true;
                else
                    fopen(obj.ViSessn);
                    pause(obj.WAIT_PAUSE_SEC);
                    ok = strcmp(obj.ViSessn.Status, 'open');                    
                end                
            catch ex
                msgString = getReport(ex);
                warning('fopen failed:\n%s',msgString);
            end
        end
		
		function Disconnect(obj)
            % Disconnect - close connection to remote instrument.
            %
            % Synopsis
            %   obj.Disconnect()
            %
            % Description
            %    Close connection to remote-instrument (if open).
            
            if strcmp(obj.ViSessn.Status, 'open')
                stopasync(obj.ViSessn);
                flushinput(obj.ViSessn);
                flushoutput(obj.ViSessn);
                fclose(obj.ViSessn);
            end
        end
        
        function [errNb, errDesc] = QuerySysErr(obj, bSendCls)
            % QuerySysErr - Query System Error from the remote instrument
            %
            % Synopsis
            %   [errNb, [errDesc]] = obj.QuerySysErr([bSendCls])
            %
            % Description
            %   Query the last system error from the remote instrument,
            %   And optionally clear the instrument's errors list.
            %
            % Inputs ([]s are optional)
            %   (bool) [bSendCls = false]  
            %           should clear the instrument's errors-list?
            %
            % Outputs ([]s are optional)
            %   (scalar) errNb     error number (zero for no error).
            %   (string) [errDesc] error description.
            
            if ~exist('bSendCls', 'var')
                bSendCls = false;
            end
            
            obj.waitTransferComplete();
            [answer, count, errmsg] = query(obj.ViSessn, 'SYST:ERR?');
            obj.waitTransferComplete();
                        
            if ~isempty(errmsg)
                error('getError() failed: %s', errmsg);
            end
            
            sep = find(answer == ',');
            if (isempty(sep) || count <= 0 || answer(count) ~= newline)
                warning('querySysErr() received invalid answer: "%s"', answer);
                flushinput(obj.ViSessn);
            end
            
            if ~isempty(sep) && isempty(errmsg)
                errNb = str2double(answer(1:sep(1) - 1));
                errmsg = answer(sep(1):end);
                if 0 ~= errNb && nargin > 1 && bSendCls
                    query(obj.ViSessn, '*CLS; *OPC?');
                end
            else
                errNb =  -1;
                if isempty(errmsg)
                    errmsg = answer;
                end               
            end
            
            if nargout > 1
                errDesc = errmsg;
            end
        end
        
        function SendCmd(obj, cmdFmt, varargin)
            % SendCmd - Send SCPI Command to instrument
            %
            % Synopsis
            %   obj.SendCmd(cmdFmt, ...)
            %
            % Description
            %   Send SCPI Command to the remote instrument.
            %
            % Inputs ([]s are optional)
            %   (string) cmdFmt      command string-format (a la printf).
            %            varargin    arguments for cmdFmt
            obj.waitTransferComplete();
            
            if nargin > 2
                cmdFmt = sprintf(cmdFmt, varargin{1:end});                
            end
            
            resp = '';
            errMsg = '';
            respLen = 0;
            
            if obj.ParanoiaLevel == 0
                fprintf(obj.ViSessn, cmdFmt);
                obj.waitTransferComplete();
            elseif obj.ParanoiaLevel == 1
                cmdFmt = strcat(cmdFmt, ';*OPC?');
                [resp, respLen, errMsg] = query(obj.ViSessn, cmdFmt);
            elseif obj.ParanoiaLevel >= 2
                cmdFmt = strcat(cmdFmt, ';:SYST:ERR?');
                [resp, respLen, errMsg] = query(obj.ViSessn, cmdFmt);
            end
            
            if (obj.ParanoiaLevel > 0 && ~isempty(errMsg))
                error('query(''%s\'') failed\n %s', cmdFmt, errMsg);
            elseif (obj.ParanoiaLevel >= 2 && respLen > 0)
                resp = deblank(resp);
                sep = find(resp == ',');
                if ~isempty(sep)
                    errNb = str2double(resp(1:sep(1) - 1));
                    if 0 ~= errNb
                        query(obj.ViSessn, '*CLS; *OPC?');
                        warning('System Error #%d after ''%s'' (%s).', ...
                            errNb, cmdFmt, resp);
                    end
                end
            end
        end
        
        
%         function [res] = SendScpi(obj, cmdFmt, varargin)
%             % SendCmd - Send SCPI Command to instrument
%             %
%             % Synopsis
%             %   obj.SendCmd(cmdFmt, ...)
%             %
%             % Description
%             %   Send SCPI Command to the remote instrument.
%             %
%             % Inputs ([]s are optional)
%             %   (string) cmdFmt      command string-format (a la printf).
%             %            varargin    arguments for cmdFmt
%             obj.waitTransferComplete();
%             
%             if nargin > 2
%                 cmdFmt = sprintf(cmdFmt, varargin{1:end});                
%             end
%             
%             res.ErrCode = 0;
%             res.RespStr = '';
%             resp = '';
%             errMsg = '';
%             respLen = 0;
%             
%             if obj.ParanoiaLevel == 0
%                 fprintf(obj.ViSessn, cmdFmt);
%                 obj.waitTransferComplete();
%             elseif obj.ParanoiaLevel == 1
%                 cmdFmt = strcat(cmdFmt, ';*OPC?');
%                 [resp, respLen, errMsg] = query(obj.ViSessn, cmdFmt);
%             elseif obj.ParanoiaLevel >= 2
%                 cmdFmt = strcat(cmdFmt, ';:SYST:ERR?');
%                 [resp, respLen, errMsg] = query(obj.ViSessn, cmdFmt);
%             end
%             
%             if (obj.ParanoiaLevel > 0 && ~isempty(errMsg))
%                 res.ErrCode = 1;
%                 error('query(''%s\'') failed\n %d %s', cmdFmt, res.ErrCode, errMsg);
%             elseif (obj.ParanoiaLevel >= 2 && respLen > 0)
%                 resp = deblank(resp);
%                 sep = find(resp == ',');
%                 if ~isempty(sep)
%                     errNb = str2double(resp(1:sep(1) - 1));
%                     if 0 ~= errNb
%                         query(obj.ViSessn, '*CLS; *OPC?');
%                         warning('System Error #%d after ''%s'' (%s).', ...
%                             errNb, cmdFmt, resp);
%                     end
%                 end
%             end
%         end
        
        
        function res = SendScpi(obj, qformat, varargin)
            % SendQuery - Send SCPI Query to instrument
            %
            % Synopsis
            %   resp = obj.SendQuery(qformat, ...)
            %
            % Description
            %   Send SCPI Query to the remote instrument,
            %   And return the instrument's response (string).
            %
            % Inputs ([]s are optional)
            %   (string) qformat     query string-format (a la printf).
            %            varargin    arguments for qformat
            %
            % Outputs ([]s are optional)
            %   (string) resp     the instrument's response.
            
            res.ErrCode = 0;
            res.RespStr = '';
            obj.waitTransferComplete();
            if nargin == 2
                if contains(qformat,'?') == 0
                    qformat = strcat(qformat, ';:SYST:ERR?');
                    [resp, respLen, errMsg] = query(obj.ViSessn, qformat);
                    respSplit = strsplit(resp,',');
                    respErr = str2double(string(respSplit(1)));
                else
                    [resp, respLen, errMsg] = query(obj.ViSessn, qformat);
                    respErr = 0;
                end
            elseif nargin > 2
                qformat = sprintf(qformat, varargin{1:end});
                [resp, respLen, errMsg] = query(obj.ViSessn, qformat);
            else
                resp = '';
                errMsg = '';
                respLen = 0;
            end
            
            
            
            if ~isempty(errMsg) || respErr~=0
                res.ErrCode = 1;
            end
            
            if res.ErrCode == 1
                error('query(''%s\'') failed\n %s', qformat, errMsg);
            end
            
            if respLen > 0
                % remove trailing blanks
                resp = deblank(resp);
                res.RespStr = resp;
            end
        end
        
        
        function res = WriteBinaryData(obj, pref, datArray, elemType, start, numItems)
            
            % SendBinaryData - Send binary data to instrument
            %
            % Synopsis
            %   obj.SendBinaryData(pref, datArray, elemType, [offset,[count]])
            %
            % Description
            %   Send array of basic-type elements to the remote instrument
            %   as binary-data with binary-data header and (optional) SCPI
            %   statement prefix (e.g. ":TRAC:DATA").
            %
            % Inputs ([]s are optional)
            %   (string) pref      SCPI statement (e.g. ":TRAC:DATA")
            %                      sent before the binary-data header.
            %   (array)  datArray  array of fixed-size elements.
            %   (string) elemType  element type name (e.g. 'uint8')
            %   (integer) [start=1] start position inside the datArray.
            %   (integer) [numItems=length(datArray)-start+1] number of
            %   elements to send from the datArray.
            
            obj.waitTransferComplete();
            
            res.ErrCode = 0;
            
            if ~exist('pref', 'var')
                pref = '';
            end            
            if ~exist('datArray', 'var')
                datArray = [];
            end            
            if ~exist('elemType', 'var')
                elemType = 'uint8';
                datArray = typecast(datArray, 'uint8');
            end
            
            if ~exist('start', 'var')
                start = 1;
            end 
            
            if ~exist('numItems', 'var')
                numItems = length(datArray) - start + 1;
            end 
            
            %numItems = length(datArray);  
            switch elemType
                case { 'int8', 'uint8' 'char' }
                    itemSz = 1;
                case { 'int16', 'uint16' }
                    itemSz = 2;
                case { 'int32', 'uint32', 'single' }
                    itemSz = 4;
                case { 'int64', 'uint64', 'double' }
                    itemSz = 8;
                otherwise
                    error('unsopported element-type ''%s''', elemType);
            end
            
            assert(itemSz >= 1 && itemSz <= obj.BINARY_CHUNK_SIZE);
            
            getChunk = @(offs, len) datArray(offs + 1 : offs + len);
            
            % make binary-data header
            szStr = sprintf('%lu', numItems * itemSz);
            pref = sprintf('*OPC?;%s#%u%s', pref, length(szStr), szStr);
            % send it (without terminating new-line!):            
            fwrite(obj.ViSessn, pref, 'char');
            obj.waitTransferComplete();
            
            % send the binary-data (in chunks):            
            
            chunkLen = fix(obj.BINARY_CHUNK_SIZE / itemSz);
            offset = 0;
            while offset < numItems
                if offset + chunkLen > numItems
                    chunkLen = numItems - offset;
                end
                dat = getChunk(start + offset - 1, chunkLen);
                fwrite(obj.ViSessn, dat, elemType);
                obj.waitTransferComplete();                
                offset = offset + chunkLen;
            end
            
            % read back the response to that *OPC? query:
            fscanf(obj.ViSessn, '%s');
            %fgets(obj.ViSessn, 2);
            
            if obj.ParanoiaLevel >= 2
                [errNb, errDesc] = obj.QuerySysErr(1);
                if 0 ~= errNb
                    res.ErrCode = 1;
                    warning('System Error #%d (%s) after ' ...
                        + 'sending ''%s ..''.', errNb, errDesc, pref);
                end
            end
        end
        
        function SendBinaryData(obj, pref, datArray, elemType, start, numItems)
            
            % SendBinaryData - Send binary data to instrument
            %
            % Synopsis
            %   obj.SendBinaryData(pref, datArray, elemType, [offset,[count]])
            %
            % Description
            %   Send array of basic-type elements to the remote instrument
            %   as binary-data with binary-data header and (optional) SCPI
            %   statement prefix (e.g. ":TRAC:DATA").
            %
            % Inputs ([]s are optional)
            %   (string) pref      SCPI statement (e.g. ":TRAC:DATA")
            %                      sent before the binary-data header.
            %   (array)  datArray  array of fixed-size elements.
            %   (string) elemType  element type name (e.g. 'uint8')
            %   (integer) [start=1] start position inside the datArray.
            %   (integer) [numItems=length(datArray)-start+1] number of
            %   elements to send from the datArray.
            
            obj.waitTransferComplete();
            
            if ~exist('pref', 'var')
                pref = '';
            end            
            if ~exist('datArray', 'var')
                datArray = [];
            end            
            if ~exist('elemType', 'var')
                elemType = 'uint8';
                datArray = typecast(datArray, 'uint8');
            end
            
            if ~exist('start', 'var')
                start = 1;
            end 
            
            if ~exist('numItems', 'var')
                numItems = length(datArray) - start + 1;
            end 
            
            %numItems = length(datArray);  
            switch elemType
                case { 'int8', 'uint8' 'char' }
                    itemSz = 1;
                case { 'int16', 'uint16' }
                    itemSz = 2;
                case { 'int32', 'uint32', 'single' }
                    itemSz = 4;
                case { 'int64', 'uint64', 'double' }
                    itemSz = 8;
                otherwise
                    error('unsopported element-type ''%s''', elemType);
            end
            
            assert(itemSz >= 1 && itemSz <= obj.BINARY_CHUNK_SIZE);
            
            getChunk = @(offs, len) datArray(offs + 1 : offs + len);
            
            % make binary-data header
            szStr = sprintf('%lu', numItems * itemSz);
            pref = sprintf('*OPC?;%s#%u%s', pref, length(szStr), szStr);
            % send it (without terminating new-line!):            
            fwrite(obj.ViSessn, pref, 'char');
            obj.waitTransferComplete();
            
            % send the binary-data (in chunks):            
            
            chunkLen = fix(obj.BINARY_CHUNK_SIZE / itemSz);
            offset = 0;
            while offset < numItems
                if offset + chunkLen > numItems
                    chunkLen = numItems - offset;
                end
                dat = getChunk(start + offset - 1, chunkLen);
                fwrite(obj.ViSessn, dat, elemType);
                obj.waitTransferComplete();                
                offset = offset + chunkLen;
            end
            
            % read back the response to that *OPC? query:
            fscanf(obj.ViSessn, '%s');
            %fgets(obj.ViSessn, 2);
            
            if obj.ParanoiaLevel >= 2
                [errNb, errDesc] = obj.QuerySysErr(1);
                if 0 ~= errNb
                    warning('System Error #%d (%s) after ' ...
                        + 'sending ''%s ..''.', errNb, errDesc, pref);
                end
            end
        end
        
        function datArray = ReadBinaryData(obj, pref, elemType)
            
            % ReadBinaryData - Read binary data from instrument
            %
            % Synopsis
            %   datArray = obj.ReadBinaryData(pref, elemType)
            %
            % Description
            %   Read array of basic-type elements from the instrument.
            %
            % Inputs ([]s are optional)
            %   (string) pref      SCPI statement (e.g. ":TRAC:DATA")
            %                      sent before the binary-data header.
            %   (string) elemType  element type name (e.g. 'uint8')
            %
            % Outputs ([]s are optional)
            %   (array)  datArray  array of fixed-size elements.
            
            obj.waitTransferComplete();
            
            if ~exist('pref', 'var')
                pref = '';
            end            
            
            switch elemType
                case { 'int8', 'uint8' 'char' }
                    itemSz = 1;
                case { 'int16', 'uint16' }
                    itemSz = 2;
                case { 'int32', 'uint32', 'single' }
                    itemSz = 4;
                case { 'int64', 'uint64', 'double' }
                    itemSz = 8;
                otherwise
                    error('unsopported element-type ''%s''', elemType);
            end
            
            assert(itemSz >= 1 && itemSz <= obj.BINARY_CHUNK_SIZE);
            
            
            % Send the prefix (if it is not empty)
            if ~isempty(pref)
                fprintf(obj.ViSessn, pref);
            end
            obj.waitTransferComplete();
            
            % Read binary header
            while true
                ch = fread(obj.ViSessn, 1, 'char');
                if ch == '#'
                    break
                end
            end
            
            % Read the first digit
            ch = fread(obj.ViSessn, 1, 'char');
            assert ('0' < ch && ch <= '9');
            
            ndigits = ch - '0';
            %fprintf('ReadBinaryData: ndigits = %d\n', ndigits);
            
            sizestr = fread(obj.ViSessn, ndigits, 'char');
            numbytes = 0;
            for n = 1:ndigits
                ch = sizestr(n, 1);
                numbytes = numbytes * 10 + (ch - '0');
            end
            
            %fprintf('ReadBinaryData: numbytes = %d\n', numbytes);
            
            datLen = ceil(numbytes / itemSz);
            assert(datLen * itemSz == numbytes);
            datArray = zeros(1, datLen, elemType);
            
            chunkLen = fix(obj.BINARY_CHUNK_SIZE / itemSz);
            
            %fprintf('ReadBinaryData: datLen=%d, chunkLen=%d\n', datLen, chunkLen);
            
            % send the binary-data (in chunks):            
            offset = 0;
            
            while offset < datLen
                if datLen - offset < chunkLen
                    chunkLen = datLen - offset;
                end
                datArray(offset + 1 : offset + chunkLen) = fread(obj.ViSessn, chunkLen, elemType);
                %obj.waitTransferComplete();                
                offset = offset + chunkLen;
            end
            
            % read the terminating newline character
            ch = fread(obj.ViSessn, 1, 'char');
            assert(ch == newline);
        end
        
        
        
    end % public methods
    
    methods (Access = private) % private methods
        
        function waitTransferComplete(obj)
            % waitTransferComplete - wait till transfer status is 'idle'
            while ~strcmp(obj.ViSessn.TransferStatus,'idle')
                pause(obj.WAIT_PAUSE_SEC);
            end
        end
    end % private methods
    
end

