function start_fig (num,dim)

fh=figure(num);clf;set(fh, 'color', 'white');
set(num,'Units', 'pixels', ...
'Position', [100 100 300 * dim(1) 300*dim(2)]);clf;

   gap=[0.05 0.05];
marg_h=0.21;
marg_v=0.14;
subtightplot(1,1,1,gap,marg_h,marg_v);
end