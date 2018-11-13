% This function is to plot needed variables

% Define plotting supports
if(strcmp(flagslice,'zonal'))
    xaxis = 0:dy:my*ly*dy;
    yaxis = 0:dz:mz*lz*dz;
    xlab = 'Horizontal distance (y) in m';
    ylab = 'Vertical distance (z) in m';
end

if(strcmp(flagslice,'meridional'))
    xaxis = 0:dx:mx*lx*dx;
    yaxis = 0:dz:mz*lz*dz;
    xlab = 'Horizontal distance (x) in m';
    ylab = 'Vertical distance (z) in m';
end

if(strcmp(flagslice,'horizontal'))
    xaxis = 0:dy:my*ly*dy;
    yaxis = 0:dx:mx*lx*dx;
    xlab = 'Horizontal distance (y) in m';
    ylab = 'Horizontal distance (x) in m';
end

if(strcmp(flagslice,'horizontalslice'))
    xaxis = 0:dy:my*ly*dy;
    yaxis = 0:dx:mx*lx*dx;
    xlab = 'Horizontal distance (y) in m';
    ylab = 'Horizontal distance (x) in m';
end

% Define what to plot

if any(strcmp(flagpar,'u'))
titlevar = ['Horizonal Fluid Velocity u_x (m/s) at time=', num2str(t),' s'];
outvar = datafullset(:,:,2)./datafullset(:,:,1);
end 

if any(strcmp(flagpar,'us'))
titlevar = ['Scaled Horizonal Fluid Velocity u_x (m/s) at time=', num2str(t),' s'];
w = datafullset(:,:,2)./datafullset(:,:,1);
for ii=1:mx*lx
for jj=1:mz*lz
outvar(ii,jj)=w(ii,jj)*scale(jj);
end
end
end

if any(strcmp(flagpar,'v'))
titlevar = ['Horizonal Fluid Velocity v_x (m/s) at time=', num2str(t),' s'];
outvar = datafullset(:,:,3)./datafullset(:,:,1);
end 

if any(strcmp(flagpar,'vs'))
titlevar = ['Scaled Horizonal Fluid Velocity v_x (m/s) at time=', num2str(t),' s'];
w = datafullset(:,:,3)./datafullset(:,:,1);
for ii=1:mx*lx
for jj=1:mz*lz
outvar(ii,jj)=w(ii,jj)*scale(jj);
end
end
end

if any(strcmp(flagpar,'w'))
titlevar = ['Vertical Fluid Velocity u_z (m/s) at time=', num2str(t),' s'];
outvar = datafullset(:,:,4)./datafullset(:,:,1);
end 

if any(strcmp(flagpar,'ws'))
titlevar = ['Scaled Vertical Fluid Velocity u_z (m/s) at time=', num2str(t),' s'];
w = datafullset(:,:,4)./datafullset(:,:,1);
for ii=1:mx*lx
for jj=1:mz*lz
outvar(ii,jj)=w(ii,jj)*scale(jj);
end
end
end

if any(strcmp(flagpar,'rhop'))
outvar=datafullset(:,:,1)-rhoa;
titlevar = ['Density Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar,'rhorp'))
outvar=100.*(datafullset(:,:,1)-rhoa)./rhoa;
titlevar = ['Density Relative Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar,'doxp'))
dox=6.022d23.*datafullset(:,:,5).*datafullset(:,:,1).*1e-3.*(1/16);
outvar=dox-doxa;
titlevar = ['Atomic Oxygen Perturbations at time=', num2str(t),' s'];
end
    
if any(strcmp(flagpar,'dnit2p'))
dnit2=6.022d23.*(1-datafullset(:,:,5)-datafullset(:,6)).*datafullset(:,1).*1e-3.*(1/28);
outvar=dnit2-dnit2a;
titlevar = ['Nitrogen Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar,'dox2p'))
dox2=6.022d23.*datafullset(:,6).*datafullset(:,1).*1e-3.*(1/32);
outvar=dox2-dox2a;
titlevar = ['Oxygen Perturbations at time=', num2str(t),' s'];
end

%    Might be useful in future 

%    % [H] density
%    dhyd=data(:,9);
%    dhyd=reshape(dhyd,mz,my,mx);
%    dhyd=permute(dhyd,[3 2 1]).*(dox+dox2+dnit2);
%    if Frame==0 dhyd0=dhyd; end;
%    dhydp=(dhyd-dhyd0);
%    
%    % [O3] density
%    dox3=data(:,10);
%    dox3=reshape(dox3,mz,my,mx);
%    dox3=permute(dox3,[3 2 1]).*(dox+dox2+dnit2);
%    if Frame==0 dox30=dox3; end;
%    dox3p=(dox3-dox30);


    if any(strcmp(flagpar,'temp'))
   rho2=datafullset(:,:,1);
   energy=datafullset(:,:,5);
   momnt=datafullset(:,:,2:4);
   momnt2=momnt.*momnt;
   kinetic=0.5*sum(momnt2,3)./rho2;
   T=(gammam-1).*((energy-kinetic)./(rho2.*Rm));
   if Frame==0 T0=T; end;
   outvar=T-T0;
   titlevar = ['Temperature pert at time=', num2str(t),' s'];
    end

        figure('pos',[500 500 600 600])
        imagesc(xaxis,yaxis,outvar')
        axis xy
        xlabel(xlab,'FontSize',14);
        ylabel(ylab,'FontSize',14);
        title(titlevar,'FontSize',14);
        colorbar 
         
                
        datafullset = [];