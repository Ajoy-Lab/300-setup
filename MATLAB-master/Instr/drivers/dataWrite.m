function retMsg = dataWrite (ip, data)
    data=uint8(data);
    t = tcpclient(ip, 65432);
    t.Timeout = 60;
    write(t, data);
    retMsg = data;
    clear t;
end