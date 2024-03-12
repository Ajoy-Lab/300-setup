function Pines_write(n)
assert(ischar(n), "n has to be integer character");
u2 = udp('192.168.1.2', 'RemotePort', 2021, 'LocalPort', 1901);
u2.EnablePortSharing = 'on';
try
    fopen(u2);
catch
    fclose(u2);
end
fwrite(u2, n);
fclose(u2);