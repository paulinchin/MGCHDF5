%%
close all
clear all
clc
it=1;
listofproc = [];

                    %-------------- Setup --------------%
                    
% Set 1 to plots figures
videooutput = 1;
% Set 1 to collect slides for creating a video
figuresoutput = 1;

% Set what input to expect.
%ouputtype = 'v'; % Vertical slice (out3ver)
%ouputtype = 'h'; % Horizontal slice (out3hor)
ouputtype = 'all'; % Full 3D domain (out3)
%ouputtype = 'air'; % Full 3D domain (out3)

% Uncomment one of the flags for needed slice:
%flagslice = 'zonal'; % constant x slice
%flagslice = 'meridional'; % constant y slice
%flagslice = 'horizontal'; % constant z slice
%flagslice = 'horizontalslice'; % constant z from horizontal sliced output
flagslice = 'full'; % constant y slice
%flagslice = 'air'; % constant y slice

geminioutput = 1;

% Set distance of slice from 0 in meters in needed direction
slicekm = 85000;

% Uncomment flags of variable to output
% Use proposed values or set index (number) of meqn (e.g. 1 - rho)
% If other variables are needed - correct plotting.m
% If scaled velocities are needed - add s to the end of variable (e.g. 'ws')
% Currently among possible outputs:
% {'u','v','w','rhop','rhorp','doxp','dnit2p','dox2p','temp'};
flagpar = {'w'};

% Do initialization
initialization

MaxFrames = 900; % Number of Frames
Frame=0; % Frame to start

% Go from frame to frame
slicing

    if(videooutput)
        % tag filename. For outstanding video quality - use ('video.avi','Uncompressed AVI')
%         v = VideoWriter('name.mp4','MPEG-4');
%         v.FrameRate = 25;
%         open(v);
%         writeVideo(v,mov(1:end));
%         close(v);
    end