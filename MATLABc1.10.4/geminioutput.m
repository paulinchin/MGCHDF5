% Subroutine to save data from MAGIC to GEMINI

velx = datafullset(:,:,2)./datafullset(:,:,1);
velz = datafullset(:,:,4)./datafullset(:,:,1);
dox=6.022d23.*datafullset(:,:,5).*datafullset(:,:,1).*1e-3.*(1/16);
dnit2=6.022d23.*(1-datafullset(:,:,5)-datafullset(:,6)).*datafullset(:,1).*1e-3.*(1/28);
dox2=6.022d23.*datafullset(:,6).*datafullset(:,1).*1e-3.*(1/32);
rho2=datafullset(:,:,1);
energy=datafullset(:,:,5);
momnt=datafullset(:,:,2:4);
momnt2=momnt.*momnt;
kinetic=0.5*sum(momnt2,3)./rho2;
T=(gammam-1).*((energy-kinetic)./(rho2.*Rm));

if (Frame == 0)
velxs = zeros(900,mx,mz);
velzs = zeros(900,mx,mz);
doxss = zeros(900,mx,mz);
dnit2s = zeros(900,mx,mz);
dox2s = zeros(900,mx,mz);
temps =zeros(900,mx,mz);
    
    velx0=velx;
    velz0=velz;
    temp0=velx;
    dox20=velx;
    dnit20=dnit2;
    dox0=dox; 
    T0 = T;
end

velxs(Frame+1,:,:)=velx;
velzs(Frame+1,:,:)=velz;
doxss(Frame+1,:,:)=dox-dox0;
dnit2s(Frame+1,:,:)=dnit2-dnit20;
dox2s(Frame+1,:,:)=dox2-dox20;
temps(Frame+1,:,:) = T-T0;
 
datafullset=[];