%% This script is to find needed threads to output from MAGIC for vertical slices
clear all
clc

% Set the values above as you use in clawez.data:
dx = 2000;
dy = 2000;
dz = 200;
lx = 10;
ly = 24;
lz = 1;
mx = 20;
my = 10;
mz = 1350;

slice = 240000; % y direction

% Set slice direction
flagslice = 'y'; % Meridional
%flagsice = 'x'; % Zonal

id = 0; % do not change this
f = 0; % do not change this

fprintf('Output these threads: \n')

if (flagslice=='y')
for ii=0:1:lx-1
for jj=0:1:ly-1
xlower = ii*mx*dx;
ylower = jj*my*dy;
xhigher = xlower + (mx)*dx;
yhigher = ylower + (my)*dy;
if ((slice>=ylower) && (slice<yhigher))
f = f+1;
id
end
id = id+1;
end
end
fprintf('Total number of cores containing this slice is: %d \n',f)
end

id = 0; % do not change this
f = 0; % do not change this

if (flagslice=='x')
for ii=0:1:lx-1
for jj=0:1:ly-1
xlower = ii*mx*dx;
ylower = jj*my*dy;
xhigher = xlower + (mx)*dx;
yhigher = ylower + (my)*dy;
if ((slice>=xlower) && (slice<xhigher))
f = f+1;
id
end
id = id+1;
end
end
fprintf('Total number of cores containing this slice is: %d \n',f)
end
