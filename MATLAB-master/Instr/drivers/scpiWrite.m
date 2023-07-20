function retMsg = scpiWrite (ip, command)
    t = tcpip(ip, 65432);
    t.Timeout = 60;
    fopen(t);
    fprintf(t, command);
    data = fscanf(t, '%s\n');
    clear t;
    retMsg = data;
end