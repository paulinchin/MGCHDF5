% set plot limits
xmin = -360;
xmax = 360;
ymin = -360;
ymax = 360;

% first pass at matlab 3d for graphics... plots isosurface.
% this needs to be improved  

 % read time and number of grids:
 n1 = Frame+1000;
 fname = ['fort.q',num2str(n1),'.hdf'];
 fname(7) = '0';
  
 if ~exist(fname) 
     disp(' ');
     disp(['Frame ',num2str(Frame),' does not exist ***']);
%     if Frame==0
%	% no initial data to plot, go on to Frame 1:
%	Frame = 1;
%	fname(length(fname))='1';
%	end
     end

if exist(fname)


% Get some data by reading first grid info 
 SD_id=hdfsd('start',fname,'read');
 [nsds,nattr,status]=hdfsd('fileinfo',SD_id); 
 sds_index=0;
 sds_id=hdfsd('select',SD_id,sds_index);
 [data,status]=hdfsd('readdata',sds_id,0,1,17); 
 status = hdfsd('endaccess',sds_id);
 status = hdfsd('end',SD_id);
 t = data(3);
 meqn = data(4);
 ngrids=nsds/(meqn+1); 

 disp(' ')
 disp(['Frame ',num2str(Frame),' at time t = ',num2str(t)]);


 
 qmin = 1e6;
 qmax = -1e6;
 ncells = [];

 if exist('beforeframe')==2
    beforeframe  % make an m-file with this name for any other commands you
                % want executed before drawing each frame, for example
                 % if you want to use axes to specify exactly where the
                 % plot will be in the window, aspect ratio, etc.
    end 


 %=============================================
 % MAIN LOOP ON GRIDS FOR THIS FRAME:
 %=============================================
 
 % Reopen the file
 SD_id=hdfsd('start',fname,'read');
 % Initialize data set counter 
 sds_index=0; 

 for ng = 1:ngrids

   % read parameters for this grid:      
   sds_id=hdfsd('select',SD_id,sds_index);
   [data,status]=hdfsd('readdata',sds_id,0,1,17);      
   status = hdfsd('endaccess',sds_id);  
   gridno=data(1); nDim=data(2); level=data(5)+1; 
   mx = data(6); my = data(7); mz=data(8);
   xlow = data(10); ylow = data(11); zlow = data(12);
   xhigh = data(14); yhigh = data(15); zhigh = data(16);
   dx = (xhigh-xlow)/mx; dy = (yhigh-ylow)/my; dz = (zhigh-zlow)/mz;

   data=zeros([mx*my*mz,meqn]);
   start=zeros([nDim 1]); stride=ones([nDim 1]); 
   edge=stride; edge(1)=mx; edge(2)=my; edge(3)=mz;
   for nq=1:meqn     
     sds_index = sds_index+1;
     sds_id=hdfsd('select',SD_id,sds_index);  
     [qcomp,status]=hdfsd('readdata',sds_id,start,stride,edge);
     status = hdfsd('endaccess',sds_id);
     data(:,nq)=reshape(qcomp,[mx*my*mz 1]);
   end      
   % Get ready for next grid   
   sds_index = sds_index+1;
   
   if Frame==0,
      if msis==1,
          pf=load('profile.data');
          whos pf
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
      else
          g=9.8;
          gamma=1.4*ones(1,mz)'; R=287*ones(1,mz)';
          rho=reshape(data(:,1),mx,my,mz);
          P=reshape(feval('pressure',data),mx,my,mz);
          rho=squeeze(rho(1,1,:));
          P=squeeze(P(1,1,:));
          T=P./(rho.*R);
          P0=max(P);
          theta=T.*(P0./P).^((gamma-1)./gamma);
          N2=g/dy.*diff(log(theta));
          %scale=(g*(gamma-1))./(gamma.*N2);
          scale=(rho./rho(1)).^(1/2);
      end;
   end;

   % xvelocity
   u=data(:,2)./data(:,1);
   u=reshape(u,mz,my,mx);
   u=permute(u,[3 2 1]);
   % matrix scaling factor, gamma, Rm
   for i=1:size(u,1), 
       scalem(i,:,:)=repmat(scale,1,size(u,2))'; 
       gammam(i,:,:)=repmat(gamma,1,size(u,2))';
       Rm(i,:,:)=repmat(R,1,size(u,2))';
   end;
   us=u.*scalem;
   
   % yvelocity
   v=data(:,3)./data(:,1);
   v=reshape(v,mz,my,mx);
   v=permute(v,[3 2 1]);
   vs=v.*scalem;

   % zvelocity
   w=data(:,4)./data(:,1);
   w=reshape(w,mz,my,mx);
   w=permute(w,[3 2 1]);
   ws=w.*scalem;

   % density
   rho=data(:,1);
   rho=reshape(rho,mz,my,mx);
   rho=permute(rho,[3 2 1]);
   if Frame==0 rho0=rho; end;
   rhop=(rho-rho0);
   rhorp=100.*(rho-rho0)./rho0;
   
   % [O] density
   dox=6.022d23.*data(:,6).*data(:,1).*1e-3.*(1/16);
   dox=reshape(dox,mz,my,mx);
   dox=permute(dox,[3 2 1]);
   if Frame==0 dox0=dox; end;
   doxp=(dox-dox0);
   
   % [N2] density
   dnit2=6.022d23.*(1-data(:,6)-data(:,7)).*data(:,1).*1e-3.*(1/28);
   dnit2=reshape(dnit2,mz,my,mx);
   dnit2=permute(dnit2,[3 2 1]);
   if Frame==0 dnit20=dnit2; end;
   dnit2p=(dnit2-dnit20);
   
   % [O2] density
   dox2=6.022d23.*data(:,7).*data(:,1).*1e-3.*(1/32);
   dox2=reshape(dox2,mz,my,mx);
   dox2=permute(dox2,[3 2 1]);
   if Frame==0 dox20=dox2; end;
   dox2p=(dox2-dox20);
   
   % [H] density
   dhyd=data(:,9);
   dhyd=reshape(dhyd,mz,my,mx);
   dhyd=permute(dhyd,[3 2 1]).*(dox+dox2+dnit2);
   if Frame==0 dhyd0=dhyd; end;
   dhydp=(dhyd-dhyd0);
   
   % [O3] density
   dox3=data(:,10);
   dox3=reshape(dox3,mz,my,mx);
   dox3=permute(dox3,[3 2 1]).*(dox+dox2+dnit2);
   if Frame==0 dox30=dox3; end;
   dox3p=(dox3-dox30);
   
   % temperature
   rho2=data(:,1);
   energy=data(:,5);
   momnt=data(:,2:4);
   energy=reshape(energy,mz,my,mx);
   momnt2=momnt.*momnt;
   kinetic=0.5*sum(momnt2,2)./rho2;
   kinetic=reshape(kinetic,mz,my,mx);
   kinetic=permute(kinetic,[3 2 1]);
   energy=permute(energy,[3 2 1]);
   T=(gammam-1).*((energy-kinetic)./(rho.*Rm));
   if Frame==0 T0=T; end;
   Tp=T-T0;
   
   % -----------------------------------------------
   % plot commands go here.... 
   if Frame>0
    x=[xlow:dx:xhigh];
    y=[ylow:dy:yhigh];
    z=[zlow:dz:zhigh];
    % set slice positions
    mpxz=uint8(1);
    mpyz=uint8(mx);
    mpxy=uint8(mz./2);
    % plot x-velocities
    figure(10);
    pdata1=squeeze(u(:,mpxz,:))';
    pdata2=squeeze(u(mpyz,:,:))';
    pdata3=squeeze(u(:,:,mpxy))';
    subplot(1,2,1);
    imagesc(x./1000,z./1000,pdata1);
    axis xy;
    axis image;
    axis([xmin xmax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    caxis([-max(max(abs(pdata1))) max(max(abs(pdata1)))]);
    colorbar;
    title('x-velocity (m/s), xz');
    xlabel('Streamwise Distance (km)');
    ylabel('Altitude (km)');
    subplot(1,2,2);
    imagesc(y./1000,z./1000,pdata2);
    axis xy;
    axis image;
    axis([ymin ymax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    caxis([-max(max(abs(pdata2))) max(max(abs(pdata2)))]); 
    colorbar;
    title('x-velocity (m/s), yz');
    xlabel('Spanwise Distance (km)');
    ylabel('Altitude (km)');
    % plot y-velocities
    figure(20);
    pdata1=squeeze(v(:,mpxz,:))';
    pdata2=squeeze(v(mpyz,:,:))';
    pdata3=squeeze(v(:,:,mpxy))';
    subplot(1,2,1);
    imagesc(x./1000,z./1000,pdata1);
    axis xy;
    axis image;
    axis([xmin xmax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    caxis([-max(max(abs(pdata1))) max(max(abs(pdata1)))]);
    colorbar;
    title('y-velocity (m/s), xz');
    xlabel('Streamwise Distance (km)');
    ylabel('Altitude (km)');
    subplot(1,2,2);
    imagesc(y./1000,z./1000,pdata2);
    axis xy;
    axis image;
    axis([ymin ymax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    caxis([-max(max(abs(pdata2))) max(max(abs(pdata2)))]);
    colorbar;
    title('y-velocity (m/s), yz');
    xlabel('Spanwise Distance (km)');
    % plot z-velocities
    figure(30);
    pdata1=squeeze(w(:,mpxz,:))';
    pdata2=squeeze(w(mpyz,:,:))';
    pdata3=squeeze(w(:,:,mpxy))';
    %if Frame==1  
    %    wstore=pdata1;
    %else
    %    wstore=[wstore pdata1];
    %end;
    subplot(1,2,1);
    imagesc(x./1000,z./1000,pdata1);
    axis xy;
    axis image;
    axis([xmin xmax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    caxis([-max(max(abs(pdata1))) max(max(abs(pdata1)))]);
    colorbar;
    title('z-velocity (m/s), xz');
    xlabel('Streamwise Distance (km)');
    ylabel('Altitude (km)');
    subplot(1,2,2);
    imagesc(y./1000,z./1000,pdata2);
    axis xy;
    axis image;
    axis([ymin ymax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    caxis([-max(max(abs(pdata2))) max(max(abs(pdata2)))]); 
    colorbar;
    title('z-velocity (m/s), yz');
    xlabel('Spanwise Distance (km)');
    ylabel('Altitude (km)');
    % plot vorticity
    if my==1,
        [X,Z]=meshgrid(x(2:end),z(2:end));
        [CURLY, CAV]=curl(X,Z,squeeze(u)',squeeze(w)');
        CURLY=CURLY';
        vorty1=abs(CURLY);
        DIV=divergence(X,Z,squeeze(u)',squeeze(w)');
        div=DIV';
    else
        [X,Y,Z]=meshgrid(x(2:end),y(2:end),z(2:end));
        [CURLX, CURLY, CURLZ, CAV]=curl(X,Y,Z,permute(u,[2,1,3]),permute(v,[2,1,3]),permute(w,[2,1,3]));
        div=divergence(X,Y,Z,permute(u,[2,1,3]),permute(v,[2,1,3]),permute(w,[2,1,3]));
        CURLX=permute(CURLX,[2,1,3]);
        CURLY=permute(CURLY,[2,1,3]);
        CURLZ=permute(CURLZ,[2,1,3]);
        div=permute(div,[2,1,3]);
        vorty2=abs(squeeze(CURLY(mpyz,:,:)));
        vorty1=abs(squeeze(CURLY(:,mpxz,:)));
        div=squeeze(div(:,mpxz,:));
    end;
    fvort=figure(40);
    subplot(1,2,1);
    imagesc(x./1000,z./1000,vorty1');
    axis xy;
    axis image;
    axis([xmin xmax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    colormap(flipud(colormap));
    caxis([0 0.25]);
    %colorbar;
    title('|\psi_y| (1/s)');
    xlabel('Streamwise (km)');
    ylabel('Altitude (km)');
    if my>1,
        subplot(1,2,2);
        imagesc(y./1000,z./1000,vorty2');
        axis xy;
        axis image;
        axis([ymin ymax 70 130]);
        colormap(lbmap(1024,'RedBlue'));
        colormap(flipud(colormap));
        caxis([0 0.25]);
        colorbar;
        title('|\psi_y| (1/s)');
        xlabel('Spanwise (km)');
        ylabel('Altitude (km)');
    end;
    n1 = Frame+1000;
    pfname = ['Vort',num2str(n1),'.png'];
    pfname(5) = '0';
    print(fvort,pfname,'-dpng','-r300');
    fdiv=figure(41);
    subplot(1,2,1);
    imagesc(x./1000,z./1000,div');
    axis xy;
    axis image;
    axis([xmin xmax 70 130]);
    colormap(lbmap(1024,'RedBlue'));
    colormap(flipud(colormap));
 %   caxis([0 0.25]);
    %colorbar;
    title('Divergence (1/s)');
    xlabel('Streamwise (km)');
    ylabel('Altitude (km)');
    if my>1,
        subplot(1,2,2);
        imagesc(y./1000,z./1000,div');
        axis xy;
        axis image;
        axis([ymin ymax 70 130]);
        colormap(lbmap(1024,'RedBlue'));
        colormap(flipud(colormap));
      %  caxis([0 0.25]);
        colorbar;
        title('Divergence (1/s)');
        xlabel('Spanwise (km)');
        ylabel('Altitude (km)');
    end;
    n1 = Frame+1000;
    pfname = ['Div',num2str(n1),'.png'];
    pfname(5) = '0';
    print(fdiv,pfname,'-dpng','-r300');
    % plot densities
    fdens=figure(50);
    pdata1=squeeze(dox(:,mpxz,:))';
    pdata2=squeeze(dox3(:,mpxz,:))';
    pdata3=squeeze(dhyd(:,mpxz,:))';
    subplot(1,3,1);
    imagesc(x./1000,z./1000,pdata1);
    axis xy;
    axis image;
    axis([xmin xmax 75 115]);
    colormap(gray(1024));
    caxis([0 4E11]);
    colorbar('horiz');
    colormap(flipud(colormap));
    title('[O] (m^{-3})');
 %   xlabel('Streamwise Distance (km)');
  %  ylabel('Altitude (km)');
    subplot(1,3,2);
    imagesc(x./1000,z./1000,pdata2);
    axis xy;
    axis image;
    axis([xmin xmax 75 115]);
    colormap(gray(1024));
    caxis([0 4E8]); 
    colorbar('horiz');
    colormap(flipud(colormap));
    title('[O_3] (m^{-3})');
  %  xlabel('Spanwise Distance (km)');
  %  ylabel('Altitude (km)');
    subplot(1,3,3);
    imagesc(x./1000,z./1000,pdata3);
    axis xy;
    axis image;
    axis([xmin xmax 75 115]);
    colormap(gray(1024));
    caxis([0 1E8]); 
    colorbar('horiz');
    colormap(flipud(colormap));
    title('[H] (m^{-3})');
    n1 = Frame+1000;
    pfname = ['Dens',num2str(n1),'.png'];
    pfname(5) = '0';    
    print(fdens,pfname,'-dpng','-r300');
   % xlabel('Spanwise Distance (km)');
   % ylabel('Altitude (km)');
   end;

   % -----------------------------------------------

   %query     % Uncomment this line to pause after plotting each subgrid
	      % Useful if you want to examine data on some subgrid

   %=============================================
   end % loop on ng (plot commands for each grid)
   %=============================================

 % done with this frame
 %status = fclose(fid);

 if exist('afterframe')==2  
    afterframe  % make an m-file with this name for any other commands you
	        % want executed at the end of drawing each frame
                % for example to change the axes, or add a curve for a
                % boundary
    end

 hold off


 status = hdfsd('end',SD_id);
 end % if exist(fname)
