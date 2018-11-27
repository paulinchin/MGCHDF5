%% This script is to find needed threads to output from MAGIC for vertical slices
clear all
clc

% Nepal
% dx = 1000;
% dy = 2000;
% dz = 250;
% lx = 18;
% ly = 20;
% lz = 1;
% mx = 990;
% my = 300;
% mz = 1800;

% Tohoku Acoustics
dx = 1000;
dy = 1000;
dz = 250;
lx = 15;
ly = 24;
lz = 1;
mx = 1200;
my = 1200;
mz = 1000;

% Palu
% dx = 1000;
% dy = 1000;
% dz = 1000;
% lx = 15;
% ly = 24;
% lz = 1;
% mx = 600;
% my = 600;
% mz = 1000;

slice = 600000; % y direction

% Set slice direction
flagslice = 'y'; % Meridional
%flagslice = 'x'; % Zonal

id = 0; % do not change this
f = 0; % do not change this

fprintf('Output these threads: \n')

if (flagslice=='y')
for ii=0:1:lx-1
for jj=0:1:ly-1
xlower = ii*(mx/lx)*dx;
ylower = jj*(my/ly)*dy;
xhigher = xlower + (mx/lx)*dx;
yhigher = ylower + (my/ly)*dy;

if ((slice>=ylower) && (slice<yhigher))
f = f+1;
id+1
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
xlower = ii*(mx/lx)*dx;
ylower = jj*(my/ly)*dy;
xhigher = xlower + (mx/lx)*dx;
yhigher = ylower + (my/ly)*dy;
if ((slice>=xlower) && (slice<xhigher))
f = f+1;
id+1
end
id = id+1;
end
end
fprintf('Total number of cores containing this slice is: %d \n',f)
end

fprintf('Do not forget to output id=1 (MASTER): %d \n',f)