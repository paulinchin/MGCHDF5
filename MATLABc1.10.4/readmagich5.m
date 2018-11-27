%%
close all
clear all
clc
it=1;
listofproc = [];
movieoutput = [];
                    %-------------- Setup --------------%
                    
% Plot figures (and saves every frame to png file)
figuresoutput = 1;
% If we create video
videooutput = 1;

% What kind of output you did?
%ouputtype = 'all'; % Full 3D domain (out3)
ouputtype = 'v'; % Vertical slice (out3ver)
%ouputtype = 'h'; % Horizontal slice (out3hor)
%ouputtype = 'air'; % Full 3D domain for airglow purposes (out3)

% What do you want to retrieve from your output?
%flagslice = 'zonal'; % constant x slice
flagslice = 'meridional'; % constant y slice
%flagslice = 'horizontal'; % constant z slice from outputtype='all'
%flagslice = 'horizontalslice'; % constant z slice from outputtype='h'
%flagslice = 'full'; % % input full 4D matrix from outputtype='all'
%flagslice = 'air'; % constant y slice

% Do you need vertical slice output to GEMINI?
geminioutput = 0;

% Set distance of slice from 0 in meters in needed direction
slicekm = 600000;

% Uncomment flags of variable to output
% Use proposed values or set index (number) of meqn (e.g. 1 - rho)
% If other variables are needed - correct plotting.m
% If scaled velocities are needed - add s to the end of variable (e.g. 'ws')
% Currently among possible outputs:
% {'u','v','w','rhop','rhorp','doxp','dnit2p','dox2p','dhydp','dox3p',temp'};
flagpar = {'u','v','w'};

MaxFrames = 900; % Number of Frames
Frame=0; % Frame to start

% Go from frame to frame
slicing

if(videooutput)
    %For outstanding video quality - use ('video.avi','Uncompressed AVI')
    v = VideoWriter('TohTh.mp4','MPEG-4');
    v.FrameRate = 25;
    open(v);
    writeVideo(v,movieoutput(1:end));
    close(v);
end