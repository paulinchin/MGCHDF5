%
% plotclaw2.m
%
% generic plotting routine for clawpack and amrclaw output in matlab
% R. J. LeVeque, 1999
%
% Various parameters are set in setplot2.m
% The default version in claw/matlab/setplot2.m can be copied to your
% directory and modified to set things up differently, or type 'k'
% at the prompt to get keyboard control and change a value.
%
%---------------------------------------------------------------------

clawdim = 2;

disp(' ')
disp('plotclaw2  plots 2d results from clawpack or amrclaw')

% set plotting parameters:
whichfile = which('setplot2');
if strcmp(whichfile,'')
    disp('*** No setplot2 file found')
else
    inp = input(['Set default plotting parameters by executing'...
		' setplot2 (y if yes)? '],'s');
    if (strcmp(inp,'y'))
       setplot2
       disp(['Executing m-script ' whichfile])
       disp(' ')
    end
end
disp(' ')

% the file setprob.m can be used to set up any necessary physical parameters
% or desired values of plotting parameters for this particular problem.

whichfile = which('setprob');
if strcmp(whichfile,'')
    %disp('*** No setprob file found')
  else
    disp(['Executing m-script ' whichfile])
    disp(' ')
    setprob
  end

%=============================================
% MAIN LOOP ON FRAMES:
%=============================================

Frame = -1;  % initialize frame counter
%Frame = 49;
while Frame <= MaxFrames

    % pause for input from user to determine if we go to next frame,
    % look at data, or skip around.  This may reset Frame counter.

if exist('NoQuery')
  if NoQuery == 1
    % set NoQuery=1 if you want the plots to be produced with no
    % query in between.  Particularly useful if you want lots of frames to
    % be printed out for an animation (put a command like makeframegif
    % in afterframe.m and set NoQuery=1)
    pause(1)
    Frame = Frame + 1;
    if Frame > MaxFrames
      break;   % break out of plotclawN after last frame
    end
 %   return
  end
else


inp = 'k';
while strcmp(inp,'k')

  inp = input(...
      ['Hit <return> for next plot, or type k, r, rr, j, i, q, or ?  '],'s');

  if strcmp(inp,'?')
    disp('  k  -- keyboard input.  Type any commands and then "return"')
    disp('  r  -- redraw current frame, without re-reading data')
    disp('  rr -- re-read current file,and redraw frame');
    disp('  j  -- jump to a particular frame')
    disp('  i  -- info about parameters and solution')
    disp('  q  -- quit')
  elseif strcmp(inp,'k')
    keyboard
  elseif strcmp(inp,'r')
    % redraw:  leave Frame counter alone
    if Frame==-1
      disp('Cannot redraw yet')
      inp = 'k';
    end
  elseif strcmp(inp,'rr')
    % redraw frame AND re-read data
    amrdata = [];
  elseif strcmp(inp,'j')
    Frame = input('Frame to jump to? ');
  elseif strcmp(inp,'i')
    if clawdim == 1
      infoplot1
      disp(' ')
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      infoframe1
    end
    if clawdim == 2
      infoplot2
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      infoframe2
    end
    if clawdim == 3
      infoplot3
      disp(' ');
      disp('hit <return> for information about this frame');
      pause
      infoframe3
    end;
    inp = 'k';
  elseif isempty(inp)
    % go to next frame
    Frame = Frame + 1;
  elseif (~strcmp(inp,'q'))
    % quit handled separately below.
    % Otherwise unrecognized input, go back and try again
    inp = 'k';
  end % if strcmp
end % while strcmp
end;

% if strcmp(inp,'q')
%   % quit now
%   break
% end
    % produce the plot:
    
    hdf_frame2

end % main loop on frames
