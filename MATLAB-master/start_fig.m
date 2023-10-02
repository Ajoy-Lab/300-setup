function start_fig (num,dim)

%set figure name
fh=figure(num);
%delete all chidren of the specified figure
clf;
%Background color is white
set(fh, 'color', 'white');
% Units are in pixels. Position [left bottom width height]
% https://www.mathworks.com/help/matlab/ref/matlab.ui.figure-properties.html

set(num,'Units', 'pixels', ...
'Position', [100 100 300 * dim(1) 300*dim(2)]);

clf;

%NOT SURE WHY WE NEED THESE LINES... ==> No subplots
gap=[0.05 0.05];
marg_h=0.21;
marg_v=0.14;
subtightplot(1,1,1,gap,marg_h,marg_v);
end