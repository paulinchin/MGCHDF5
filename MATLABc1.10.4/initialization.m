    if(ouputtype=='v')
        % Vertical slice at 1 particular longitude/latitude
       filesbasement = 'fort.qv'; 
    elseif(ouputtype=='h')
        % Horizontal slice at 1 particular altitude
       filesbasement = 'fort.qh';
    elseif(ouputtype=='all')
        % Full 3D input
       filesbasement = 'fort.q'; 
    elseif(ouputtype=='air')
        % Range of horizontal slices for airglow calculation
       filesbasement = 'fort.qa';
    end
    
    nameCur = strcat(filesbasement,num2str(Frame,'%04.f'),'.h5');
    attr = h5readatt(nameCur,'/Pid1','Parameters'); % Thread MASTER (0) always outputs its data
    ngrids_out = attr(1);
    nDim = attr(2);
    t = attr(3);
    meqn = attr(4);
    mx = attr(6);
    my = attr(7);
    mz = attr(8);
    domain = [attr(10);attr(11);attr(12);attr(14);attr(15);attr(16)];
    lx = attr(22);
    ly = attr(23);
    lz = attr(24);
    iframe = attr(25);
    dx = attr(18);
    dy = attr(19);
    dz = attr(20);
    mxp = mx/lx;
    myp = my/ly;
    mzp = mz/lz;
    
disp('---------- Simulation info -----------')
fprintf('Number of cells in the domain: %dx, %dy, %dz\n',mx,my,mz);
fprintf('Physical domain: %dm, %dm, %dm\n',mx*dx,my*dy,mz*dz);
fprintf('Physical grid size: %dm, %dm, %dm\n',dx,dy,dz);
fprintf('Number of processors used: %d (%d,%d,%d)\n',lx*ly*lz,lx,ly,lz);
%fprintf('Time step: %d s\n',dt);
disp('--------------------------------------')

%-------------- Load from file .  --------------%
datafullset = hdf5read(nameCur,'/Pid1');
disp('Reading atmospheric profile...');
pf=load('profile.data');
% height

    zlower = domain(3);
    for k=-1:1:(mz+2)
    height(k+2) = zlower + (k-0.5d0)*dz;
    end

%height=pf(:,1);
%height=height(3:mz+2);
% initial total density
rho0 = squeeze(datafullset(1,1,:,1));
% gravity
g=-3.989e14./((6.38d6+height.*1000).^2);
gmid=0.5*(g(1:end-1)+g(2:end));
% initial species densities
dox0=squeeze(6.022d23.*datafullset(1,1,:,5).*datafullset(1,1,:,1).*1e-3.*(1/16));
dnit20=squeeze(6.022d23.*(1-datafullset(1,1,:,5)-datafullset(1,1,:,6)).*datafullset(1,1,:,1).*1e-3.*(1/28));
dox20=squeeze(6.022d23.*datafullset(1,1,:,6).*datafullset(1,1,:,1).*1e-3.*(1/32));
dhyd0=squeeze(datafullset(1,1,:,9));
dox30=squeeze(datafullset(1,1,:,10));
% R divided by mean molecular mass
R=8.31d3./((1d-6).*((dox0.*16+dnit20.*28+dox20.*32)./(1d-6.*(dox0+dnit20+dox20))));
% gamma
gamma=(7.*(dox20+dnit20)+5.*(dox0))./(5.*(dox20+dnit20)+3.*(dox0));
% scale factor
scale=(rho0./rho0(1)).^(1/2);
% temperature
energy=squeeze(datafullset(1,1,:,5));
momnt=squeeze(datafullset(1,1,:,2:4));
momnt2=momnt.*momnt;
kinetic=0.5*sum(momnt2,2)./rho0;
T0=(gamma-1).*((energy-kinetic)./(rho0.*R));
              
        if(strcmp(flagslice,'meridional'))
        scale = repmat(scale',mx,1);
        gamma=repmat(gamma',mx,1);
        R=repmat(R',mx,1);    
        rho0 = repmat(rho0',mx,1);
        dox0 = repmat(dox0',mx,1);
        dnit20 = repmat(dnit20',mx,1);
        dox20 = repmat(dox20',mx,1);
        dhyd0 = repmat(dhyd0',mx,1);
        dox30 = repmat(dox30',mx,1);
        end
        
        if(strcmp(flagslice,'zonal'))
        scale = repmat(scale',my,1);
        gamma=repmat(gamma',my,1);
        R=repmat(R',my,1);    
        rho0 = repmat(rho0',my,1);
        dox0 = repmat(dox0',my,1);
        dnit20 = repmat(dnit20',my,1);
        dox20 = repmat(dox20',my,1);
        dhyd0 = repmat(dhyd0',my,1);
        dox30 = repmat(dox30',my,1);
        end