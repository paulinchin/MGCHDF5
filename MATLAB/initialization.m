try
    % This script should be improved later to read simulation parameters
attr = h5readatt('fort.q0000.h5','/Pid1','Parameters');
attr2 = h5readatt('fort.q0001.h5','/Pid1','Parameters');

gridno = attr(1);
level = attr(2);
mx = attr(3);
my = attr(4);
mz = attr(5);
dx = attr(9);
dy = attr(10);
dz = attr(11);
t = attr(12);
iframe = attr(13);
meqn = attr(14);
lx = attr(15);
ly = attr(16);
lz = attr(17);
dt = attr2(12)-attr(12);

disp('---------- Simulation info -----------')
fprintf('Number of cells in the domain: %dx, %dy, %dz\n',mx*lx,my*ly,mz*lz);
fprintf('Physical domain: %dm, %dm, %dm\n',mx*lx*dx,my*ly*dy,mz*lz*dz);
fprintf('Physical grid size: %dm, %dm, %dm\n',dx,dy,dz);
fprintf('Number of processors used: %d (%d,%d,%d)\n',lx*ly*lz,lx,ly,lz);
fprintf('Time step: %d s\n',dt);
disp('--------------------------------------')

%-------------- Load MSIS profile --------------%
        %if (~exist('pf','var'))
        disp('Reading atmospheric profile...');
        attr = h5readatt('fort.q0000.h5','/Pid1','Parameters');
        pf=load('profile.data');
          %temperature
          T=pf(:,6);
          %density
          rho=pf(:,5)*1000;
          %height
          height=pf(:,1);
          %gravity
          g=-3.989e14./((6.38d6+height.*1000).^2);
          gmid=0.5*(g(1:end-1)+g(2:end));
          %atomic o
          dox0=pf(:,2);
          %molec n
          dnit20=pf(:,3);
          %molec o2
          dox20=pf(:,4);
          %R divided by mean molecular mass
          R=8.31d3./((1d-6).*((dox0.*16+dnit20.*28+dox20.*32)./(1d-6.*(dox0+dnit20+dox20))));
          %pressure
          P=rho.*R.*T;
          P0=max(P);
          %gamma
          %gamma=1.4.*(dox20+dnit20)./(dox0+dox20+dnit20)+1.67.*(dox0)./(dox0+dox20+dnit20);
          gamma=(7.*(dox20+dnit20)+5.*(dox0))./(5.*(dox20+dnit20)+3.*(dox0));
          %potential temperature
          theta=T.*(P0./P).^((gamma-1)./gamma);
          %square of buoyancy frequency
          N2=gmid./dz.*diff(log(theta));
          %fit stuff to usable grid
          P=P(3:mz+2);
          rho=rho(3:mz+2);
          height=height(3:mz+2);
          R=R(3:mz+2);
          dox0=dox0(3:mz+2);
          dnit20=dnit20(3:mz+2);
          dox20=dox20(3:mz+2);
          gamma=gamma(3:mz+2);
          T=T(3:mz+2);
          theta=theta(3:mz+2);
          N2=.5*(N2(3:mz+2)+N2(2:mz+1));
          scale=(rho./rho(1)).^(1/2);
          rho0=rho;

        if(strcmp(flagslice,'meridional'))
        rhoa = repmat(rho0',mx*lx,1);
        doxa = repmat(dox0',mx*lx,1);
        dnit2a = repmat(dnit20',mx*lx,1);
        dox2a = repmat(dox20',mx*lx,1);
        end
        
        if(strcmp(flagslice,'zonal'))
        rhoa = repmat(rho0',my*ly,1);
        doxa = repmat(dox0',my*ly,1);
        dnit2a = repmat(dnit20',my*ly,1);
        dox2a = repmat(dox20',my*ly,1);
        end

        %end
          
catch
    fprintf('fort.q0000.h5 and fort.q0001.h5 used to define simulation parameters were not found \n');
end