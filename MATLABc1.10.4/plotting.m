% Plot needed variables
% Define plotting supports
if(strcmp(flagslice,'zonal'))
    xaxis = domain(2):dy:domain(5);
    yaxis = domain(3):dz:domain(6);
    xlab = 'Horizontal distance (y) in m';
    ylab = 'Vertical distance (z) in m';
end

if(strcmp(flagslice,'meridional'))
    xaxis = domain(1):dx:domain(4);
    yaxis = domain(3):dz:domain(6);
    xlab = 'Horizontal distance (x) in m';
    ylab = 'Vertical distance (z) in m';
end

if(strcmp(flagslice,'horizontal'))
    xaxis = domain(2):dy:domain(5);
    yaxis = domain(1):dz:domain(4);
    xlab = 'Horizontal distance (y) in m';
    ylab = 'Horizontal distance (x) in m';
end

if(strcmp(flagslice,'horizontalslice'))
    xaxis = domain(2):dy:domain(5);
    yaxis = domain(1):dz:domain(4);
    xlab = 'Horizontal distance (y) in m';
    ylab = 'Horizontal distance (x) in m';
end

nsubplots = numel(flagpar);
figure('pos',[500 500 700*nsubplots 600])     


for i=1:1:nsubplots
    
subplot(1,nsubplots,i);
        
if any(strcmp(flagpar(i),'u'))
titlevar = ['Horizonal Fluid Velocity u_x (m/s) at time=', num2str(t),' s'];
outvar = datafullset(:,:,2)./datafullset(:,:,1);
end 

if any(strcmp(flagpar(i),'us'))
titlevar = ['Scaled Horizonal Fluid Velocity u_x (m/s) at time=', num2str(t),' s'];
u = datafullset(:,:,2)./datafullset(:,:,1);
outvar=u.*scale;
end

if any(strcmp(flagpar(i),'v'))
titlevar = ['Horizonal Fluid Velocity u_y (m/s) at time=', num2str(t),' s'];
outvar = datafullset(:,:,3)./datafullset(:,:,1);
end 

if any(strcmp(flagpar(i),'vs'))
titlevar = ['Scaled Horizonal Fluid Velocity u_y (m/s) at time=', num2str(t),' s'];
v = datafullset(:,:,3)./datafullset(:,:,1);
outvar=v.*scale;
end

if any(strcmp(flagpar(i),'w'))
titlevar = ['Vertical Fluid Velocity u_z (m/s) at time=', num2str(t),' s'];
outvar = datafullset(:,:,4)./datafullset(:,:,1);
end 

if any(strcmp(flagpar(i),'ws'))
titlevar = ['Scaled Vertical Fluid Velocity u_z (m/s) at time=', num2str(t),' s'];
w = datafullset(:,:,4)./datafullset(:,:,1);
outvar=w.*scale;
end

if any(strcmp(flagpar(i),'rhop'))
outvar=datafullset(:,:,1)-rho0;
titlevar = ['Density Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar(i),'rhorp'))
outvar=100.*(datafullset(:,:,1)-rho0)./rho0;
titlevar = ['Density Relative Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar(i),'doxp'))
dox=6.022d23.*datafullset(:,:,5).*datafullset(:,:,1).*1e-3.*(1/16);
outvar=dox-dox0;
titlevar = ['Atomic Oxygen Perturbations at time=', num2str(t),' s'];
end
    
if any(strcmp(flagpar(i),'dnit2p'))
dnit2=6.022d23.*(1-datafullset(:,:,5)-datafullset(:,6)).*datafullset(:,1).*1e-3.*(1/28);
outvar=dnit2-dnit20;
titlevar = ['Nitrogen Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar(i),'dox2p'))
dox2=6.022d23.*datafullset(:,6).*datafullset(:,1).*1e-3.*(1/32);
outvar=dox2-dox20;
titlevar = ['Oxygen Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar(i),'dhydp'))
dhyd=datafullset(:,9);
outvar=dhyd-dhyd0;
titlevar = ['Hydrogen Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar(i),'dox3p'))
dox3=datafullset(:,10);
outvar=dhyd-dhyd0;
titlevar = ['Ozone Perturbations at time=', num2str(t),' s'];
end

if any(strcmp(flagpar(i),'temp'))
rho2=datafullset(:,:,1);
energy=datafullset(:,:,5);
momnt=datafullset(:,:,2:4);
momnt2=momnt.*momnt;
kinetic=0.5*sum(momnt2,3)./rho2;
T=(gammam-1).*((energy-kinetic)./(rho2.*Rm));
outvar=T-T0;
titlevar = ['Temperature pert at time=', num2str(t),' s'];
end

        imagesc(xaxis,yaxis,outvar')
        axis xy
        xlabel(xlab,'FontSize',14);
        ylabel(ylab,'FontSize',14);
        title(titlevar,'FontSize',14);
        caxis([-max(max(outvar)) max(max(outvar))]);
        colorbar   
        
end 
    n1 = Frame+1000;
    pfname = ['out',num2str(n1),'.png'];
    pfname(5) = '0';
    print(pfname,'-dpng','-r300');
    
pause

datafullset = [];