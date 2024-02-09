%  Maxwell McAllister, 10/15/2020
%  Receive script, for Sage, aka Host 2
%  For communication between PCs, replace '127.0.0.1' with
% the DNS or DHCP (?) host name or IP address of Host 1
clear;
u2 = udp('192.168.1.2', 'RemotePort', 2020, 'LocalPort', 9090);
u2.EnablePortSharing = 'on';
fopen(u2);
% Loop to repeatedly wait for messages and send replies
% Break or Ctrl+C to get out of loop
while ( true )
    % Wait for message
    while(u2.BytesAvailable == 0)
        % If no bytes in u2 buffer, wait 10ms then check again
        pause(0.01);
    end
    % Bytes no available - a message has arrived
    host1Bytes = fscanf(u2);
    % Immediately send a reply
    replyTo1 = (host1Bytes + 20);
    %fwrite(u2, replyTo1);
    % Display received message and sent reply
    %disp('Message:')
    disp( host1Bytes )
    %disp('Replied:')
    %disp( replyTo1 )
    %disp('')
end
fclose(u2);
% Please execute after this script: fclose(u2)
