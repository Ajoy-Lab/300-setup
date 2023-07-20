disp('Client side, pushing data to server');
data = 'Message from client';
t = tcpip('192.168.1.2', 30000, 'NetworkRole', 'client');
fopen(t);
fwrite(t, data);
