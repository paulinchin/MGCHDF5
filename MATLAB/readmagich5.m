%%
% Function
% 1. Get information about simulation
% 2. Do vertical slice (zonal, horizontal)
% 3. Do horizontal slice
% 4. Load all matrix
a = h5info('fort.q0043.h5');
tempp = hdf5read('fort.q0043.h5','/Pid157');
b = squeeze(tempp(:,5,:,1));
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

% Set the distance of slice from (0,0) in m:
%slicex= 102000; % x direction
slicey = 240000; % y direction

calcsliceh5(slicey,dx,dy,dz,lx,ly,lz,mx,my,mz,'y')
%%
close all
clear all
clc
%-------------- Setup --------------%
% Set 1 to collect slides for creating a video
videooutput = 1;
figuresoutput = 1;

% Uncomment one of the flag for slices
%flagslice = 'zonal'; % constant x
flagslice = 'meridional'; % constant y
%flagslice = 'horizontal'; % constant z
%flagslice = 'horizontalslice'; % constant z from horizontal sliced output

% Uncomment needed distance of slice from 0 in meters
%slicezonalkm = 102000; % x direction
slicemeridionalkm = 240000; % y direction
%slicehorizontalkm = 250000; % z direction

% Uncomment needed to output flags
% Use proposed values or set index (number) of meqn (e.g. 1 - rho)
% If scaled velocities are needed - add s to the end of variable (e.g. ws)
% Currently among possible outputs:
% flagpar = {'u','v','w','rhop','rhorp','doxp','dnit2p','dox2p'};
% If other variables are needed - correct plotting.m
flagpar = {'w'};

% Do initialization
initialization

it=1;
listofproc = [];
MaxFrames = 900; % Number of Frames
Frame=132; % Frame to start

% Go from frame to frame
slicing

    if(videooutput)
        % tag filename. For outstanding video quality - use ('video.avi','Uncompressed AVI')
%         v = VideoWriter('video2.mp4','MPEG-4');
%         v.FrameRate = 25;
%         open(v);
%         writeVideo(v,mov(1:end));
%         close(v);
    end