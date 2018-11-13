% Version working with hdf5 output routines
% Now these values are set manually. Update is soon
lx = 15;
ly = 24;
mz = 50;
mx = 1200;
my = 1200;

%Plot position (grid points)
x = 1:1:1200;
y = 1:1:1200;
minX=1; 
maxX=1200;
minY=1;
maxY=1200;
sliceX=180;
sliceY=180;
%Altitude range (km)
%maxH=120;
%minH=70;
z = 76:1:125;
maxH=50;%int32(maxH*1000/dz);
minH=1;%int32(minH*1000/dz);

%First frame initialzation
if Frame==0 
  coefficients
  pf=load('profile.data');
  height=pf(:,1);
  height = height(73:122);

dox0=6.022d23.*datafullset(:,:,:,6).*datafullset(:,:,:,1).*1e-3.*(1/16);
dnit20=6.022d23.*(1-datafullset(:,:,:,6)-datafullset(:,:,:,7)).*datafullset(:,:,:,1).*1e-3.*(1/28);
dox20=6.022d23.*datafullset(:,:,:,7).*datafullset(:,:,:,1).*1e-3.*(1/32);
dhyd0 = datafullset(:,:,:,9).*(dox0+dox20+dnit20);
dox30 = datafullset(:,:,:,10).*(dox0+dox20+dnit20);
gammam=(7.*(dox20+dnit20)+5.*(dox0))./(5.*(dox20+dnit20)+3.*(dox0));
Rm=8.31d3./((1d-6).*((dox0.*16+dnit20.*28+dox20.*32)./(1d-6.*(dox0+dnit20+dox20))));
rho2=datafullset(:,:,:,1);
energy=datafullset(:,:,:,5);
momnt=datafullset(:,:,:,2:4);
momnt2=momnt.*momnt;
kinetic=0.5*sum(momnt2,4)./rho2;
T0=(gammam-1).*((energy-kinetic)./(rho2.*Rm));

    %Plot Species
    figure(101);
    sqdox0=abs(squeeze(dox0(1,1,:)));
    sqdox20=abs(squeeze(dox20(1,1,:)));
    sqdnit20=abs(squeeze(dnit20(1,1,:)));
    sqdox30=abs(squeeze(dox30(1,1,:)));
    sqdhyd0=abs(squeeze(dhyd0(1,1,:)));
    
    semilogx(sqdox0,height,'k-',sqdox20,height,'k-.',sqdnit20,height,'k:');
    legend('O','O_2','N_2');
    xlabel('Log Density [Molecules cm^-3]');
    ylabel('Altitude [km]');
    axis([10^10 10^15 75 110]);
    title('Major Species Densities');
    
    figure(102);
    semilogx(sqdox30,height,'k-',sqdhyd0,height,'k-.');
    legend('O_{3}','H');
    xlabel('Log Density [Molecules cm^-3]');
    ylabel('Altitude [km]');
    axis([10^3 10^10 75 110]);
    title('Minor Species Densities'); 
end

dox=6.022d23.*datafullset(:,:,:,6).*datafullset(:,:,:,1).*1e-3.*(1/16);
dnit2=6.022d23.*(1-datafullset(:,:,:,6)-datafullset(:,:,:,7)).*datafullset(:,:,:,1).*1e-3.*(1/28);
dox2=6.022d23.*datafullset(:,:,:,7).*datafullset(:,:,:,1).*1e-3.*(1/32);
dhyd = datafullset(:,:,:,9).*(dox+dox2+dnit2);
dox3 = datafullset(:,:,:,10).*(dox+dox2+dnit2);
rho2=datafullset(:,:,:,1);
energy=datafullset(:,:,:,5);
momnt=datafullset(:,:,:,2:4);
momnt2=momnt.*momnt;
kinetic=0.5*sum(momnt2,4)./rho2;
T=(gammam-1).*((energy-kinetic)./(rho2.*Rm));

Tp = T-T0;

    %Calculated steadz-state Makhlouf et al. OH(v) concentrations
        rk1=(1.4d-10).*exp(-470./T);
        rk8=ones(mx,my,mz).*1d-11;
        rk9=ones(mx,my,mz).*1d-13;
        rk10=ones(mx,my,mz).*1d-14;
             
%       # Calculate [OH(v=9)]
             rP=rk1.*b(9+1).*dhyd.*dox3;
             rL=rk8.*a8(9+1).*dox+rk9.*sum(a9(9,:)).*dox2+rk10.*a10(9+1).*dnit2+ones(mx,my,mz)*(aa(9,1)+aa(9,2)+aa(9,3)+aa(9,4)+aa(9,5)+aa(9,6));
             dOH9=rP./rL;
%       # Calculate [OH(v=8)]
             rP=rk1.*b(8+1).*dhyd.*dox3+dOH9.*aa(9,1)+rk9.*(a9(9,8+1).*dOH9).*dox2+rk10.*a10(9+1).*dnit2.*dOH9;
             rL=rk8.*a8(8+1).*dox+rk9.*sum(a9(8,:)).*dox2+rk10.*a10(8+1).*dnit2+ones(mx,my,mz)*(aa(8,1)+aa(8,2)+aa(8,3)+aa(8,4)+aa(8,5)+aa(8,6));
             dOH8=rP./rL;
%       # Calculate [OH(v=7)]
             rP=rk1.*b(7+1).*dhyd.*dox3+dOH9.*aa(9,2)+dOH8.*aa(8,1)+rk9.*(a9(9,7+1).*dOH9+a9(8,7+1).*dOH8).*dox2+rk10.*a10(8+1).*dnit2.*dOH8;
             rL=rk8.*a8(7+1).*dox+rk9.*sum(a9(7,:)).*dox2+rk10.*a10(7+1).*dnit2+ones(mx,my,mz)*(aa(7,1)+aa(7,2)+aa(7,3)+aa(7,4)+aa(7,5)+aa(7,6));
             dOH7=rP./rL;
%       # Calculate [OH(v=6)]
             rP=rk1.*b(6+1).*dhyd.*dox3+dOH9.*aa(9,3)+dOH8.*aa(8,2)+dOH7.*aa(7,1)+rk9.*(a9(9,6+1).*dOH9+a9(8,6+1).*dOH8+a9(7,6+1).*dOH7).*dox2+rk10.*a10(7+1).*dnit2.*dOH7;
             rL=rk8.*a8(6+1).*dox+rk9.*sum(a9(6,:)).*dox2+rk10.*a10(6+1).*dnit2+ones(mx,my,mz)*(aa(6,1)+aa(6,2)+aa(6,3)+aa(6,4)+aa(6,5)+aa(6,6));
             dOH6=rP./rL;
%       # Calculate [OH(v=5)]
             rP=dOH9.*aa(9,4)+dOH8.*aa(8,3)+dOH7.*aa(7,2)+dOH6.*aa(6,1)+rk9.*(a9(9,5+1).*dOH9+a9(8,5+1).*dOH8+a9(7,5+1).*dOH7+a9(6,5+1).*dOH6).*dox2+rk10.*a10(6+1).*dnit2.*dOH6;             
             rL=rk8.*a8(5+1).*dox+rk9.*sum(a9(5,:)).*dox2+rk10.*a10(5+1).*dnit2+ones(mx,my,mz)*(aa(5,1)+aa(5,2)+aa(5,3)+aa(5,4)+aa(5,5));
             dOH5=rP./rL;
%       # Calculate [OH(v=4)]
             rP=dOH9.*aa(9,5)+dOH8.*aa(8,4)+dOH7.*aa(7,3)+dOH6.*aa(6,2)+dOH5.*aa(5,1)+rk9.*(a9(9,4+1).*dOH9+a9(8,4+1).*dOH8+a9(7,4+1).*dOH7+a9(6,4+1).*dOH6+a9(5,4+1).*dOH5).*dox2+rk10.*a10(5+1).*dnit2.*dOH5;
             rL=rk8.*a8(4+1).*dox+rk9.*sum(a9(4,:)).*dox2+rk10.*a10(4+1).*dnit2+ones(mx,my,mz)*(aa(4,1)+aa(4,2)+aa(4,3)+aa(4,4));
             dOH4=rP./rL;
%       # Calculate [OH(v=3)]
             rP=dOH9.*aa(9,6)+dOH8.*aa(8,5)+dOH7.*aa(7,4)+dOH6.*aa(6,3)+dOH5.*aa(5,2)+dOH4.*aa(4,1)+rk9.*(a9(9,3+1).*dOH9+a9(8,3+1).*dOH8+a9(7,3+1).*dOH7+a9(6,3+1).*dOH6+a9(5,3+1).*dOH5+a9(4,3+1).*dOH4).*dox2+rk10.*a10(4+1).*dnit2.*dOH4;
             rL=rk8.*a8(3+1).*dox+rk9.*sum(a9(3,:)).*dox2+rk10.*a10(3+1).*dnit2+ones(mx,my,mz)*(aa(3,1)+aa(3,2)+aa(3,3));
             dOH3=rP./rL;
%       # Calculate [OH(v=2)]
             rP=dOH8.*aa(8,6)+dOH7.*aa(7,5)+dOH6.*aa(6,4)+dOH5.*aa(5,3)+dOH4.*aa(4,2)+dOH3.*aa(3,1)+rk9.*(a9(9,2+1).*dOH9+a9(8,2+1).*dOH8+a9(7,2+1).*dOH7+a9(6,2+1).*dOH6+a9(5,2+1).*dOH5+a9(4,2+1).*dOH4+a9(3,2+1).*dOH3).*dox2+rk10.*a10(3+1).*dnit2.*dOH3;
             rL=rk8.*a8(2+1).*dox+rk9.*sum(a9(2,:)).*dox2+rk10.*a10(2+1).*dnit2+ones(mx,my,mz)*(aa(2,1)+aa(2,2));
             dOH2=rP./rL;
%       # Calculate [OH(v=1)]
             rP=dOH7.*aa(7,6)+dOH6.*aa(6,5)+dOH5.*aa(5,4)+dOH4.*aa(4,3)+dOH3.*aa(3,2)+dOH2.*aa(2,1)+rk9.*(a9(9,1+1).*dOH9+a9(8,1+1).*dOH8+a9(7,1+1).*dOH7+a9(6,1+1).*dOH6+a9(5,1+1).*dOH5+a9(4,1+1).*dOH4+a9(3,1+1).*dOH3+a9(2,1+1).*dOH2).*dox2+rk10.*a10(2+1).*dnit2.*dOH2;
             rL=rk8.*a8(1+1).*dox+rk9.*sum(a9(1,:)).*dox2+rk10.*a10(1+1).*dnit2+ones(mx,my,mz)*aa(1,1);
             dOH1=rP./rL;
%       # Calculate [OH(v=0)]
             rP=dOH6.*aa(6,6)+dOH5.*aa(5,5)+dOH4.*aa(4,4)+dOH3.*aa(3,3)+dOH2.*aa(2,2)+dOH1.*aa(1,1)+rk9.*(a9(9,0+1).*dOH9+a9(8,0+1).*dOH8+a9(7,0+1).*dOH7+a9(6,0+1).*dOH6+a9(5,0+1).*dOH5+a9(4,0+1).*dOH4+a9(3,0+1).*dOH3+a9(2,0+1).*dOH2+a9(1,0+1).*dOH1).*dox2+rk10.*a10(1+1).*dnit2.*dOH1;
             rL=rk8.*a8(0+1).*dox;
             dOH0=rP./rL;
     
%First frame storage / plotting
if Frame==0
    %OH(v) Bands
    dOH90=dOH9;
    dOH80=dOH8;
    dOH70=dOH7;
    dOH60=dOH6;
    dOH50=dOH5;
    dOH40=dOH4;
    dOH30=dOH3;
    dOH20=dOH2;
    dOH10=dOH1;
    dOH00=dOH0;
    figure(103);
    subplot(1,2,1);
    plot(abs(squeeze(dOH90(1,1,:))),height,abs(squeeze(dOH80(1,1,:))),height,abs(squeeze(dOH70(1,1,:))),height,abs(squeeze(dOH60(1,1,:))),height);
    legend('v=9','v=8','v=7','v=6');
    axis([10^1 10^3 75 100]);
    xlabel('Concentration [cm^-3]');
    ylabel('Altitude [km]');
    title('Hydroxyl OH(v=6-9)');
    subplot(1,2,2);
    plot(abs(squeeze(dOH50(1,1,:))),height,abs(squeeze(dOH40(1,1,:))),height,abs(squeeze(dOH30(1,1,:))),height,abs(squeeze(dOH20(1,1,:))),height,abs(squeeze(dOH10(1,1,:))),height);
    legend('v=5','v=4','v=3','v=2','v=1');
    axis([10^1 10^4 75 100]);
    xlabel('Concentration [cm^-3]');
    ylabel('Altitude [km]');
    title('Hydroxyl OH(v=1-5)');

    %OH(8,3)
    figure(104);
    sqT=squeeze(T0(1,1,:));
    k1o3=1.4E-10.*exp(-470./sqT);
    k6N2=5.7E-34.*(sqT./300).^(-2.62);
    k6O2=5.96E-34.*(sqT./300).^(-2.37);
    subplot(1,2,2);
    %visOH8o=.29.*dox0.*dox20.*(k6N2.*dnit20+k6O2.*dox20)./(260+2E-11.*dox20);
    visOH8o=.29.*sqdhyd0.*sqdox30.*k1o3./(260+2E-11.*sqdox20);
    semilogx(abs(visOH8o),height,'k-',abs(.78.*squeeze(dOH80(1,1,:))),height,'k-.');
    legend('McDade et al.','Makhlouf et al.');
    axis([10^0 10^3 75 110]);
    xlabel('Volume Emission Rate [cm^-3 s^-1]');
    ylabel('Altitude [km]');
    title('Hydroxyl (8,3) Meinel Band');
    %O1
    figure(105);
    k1=4.7E-33.*(300./sqT).^2;
    k5=4.0E-12.*exp(-865./sqT);
    A5=1.18;
    A6=1.35;
    CpO=2.1E2;
    CpO2=1.5E1;
    M=sqdox20+sqdnit20+sqdox0;
    visOIo=A5.*k1.*(sqdox0.^3).*M./((A6+k5.*sqdox20).*(CpO2.*sqdox20+CpO.*sqdox0));
    %Schubert
    sk1=0.03.*k1;
    k2=5.0e-13;
    k3=3.0e-11;
    A1=2.0e-2;
    A2=1.105;
    sk3=.2*k3;
    k6=k5;
    A5577=1.06;
    nO2c=sk1.*(sqdox0.^2).*M./(k2.*sqdox20+k3.*sqdox0+A1);
    nO1s=sk3.*sqdox0.*nO2c./(k6.*sqdox20+A2);
    visOIso=A5577.*nO1s;
    %Modified Schubert / Hickey
    sk1=0.03.*k1;
    k2=5.0e-13;
    A1=2.0e-2;
    A2=1.35;
%   sk3=2.5e-12;
    k6=k5;
    A5577=1.18;
    nO2c=sk1.*(sqdox0.^2).*M./(k2.*sqdox20+k3.*sqdox0+A1);
    nO1s=sk3.*sqdox0.*nO2c./(k6.*sqdox20+A2);
    visOImso=A5577.*nO1s;
    %Modified Hickey / Schubert
    sk1=0.03.*k1;
    k2=5.0e-13;
    k3=6e-12;
    A1=1.0e-3;
    A2=1.35;
    sk3=0.2*k3;
    k6=k5;
    A5577=1.18;
    nO2c=sk1.*(sqdox0.^2).*M./(k2.*sqdox20+k3.*sqdox0+A1);
    nO1s=sk3.*sqdox0.*nO2c./(k6.*sqdox20+A2);
    visOIhso=A5577.*nO1s;
    semilogx(abs(visOIo),height,'k-',abs(visOIso),height,'k--',abs(visOImso),height,'k-.',abs(visOIhso),height,'k:');
    legend('McDade et al.','Schubert et al.','Mod. Schubert et al.','Mod. Hickey et al.');
    axis([10^0 10^3 75 110]);
    xlabel('Volume Emission Rate (cm^-3 s^-1)');
    ylabel('Altitude (km)');
    title('OI Greenline 557.7nm');      
    %O2
    figure(106);
    k1o2a=(4.7e-33).*(300./sqT).^2;
    k2n2a=2.2e-15;
    k2o2a=4e-17;
    A1=0.079; A2=0.083;
    visO2o=(k1o2a.*A1.*(sqdox0.^2).*(sqdox20+sqdnit20).*sqdox20)./((A2+k2o2a.*sqdox20+k2n2a.*sqdnit20).*(7.5.*sqdox20+33.*sqdox0));
    semilogx(abs(visO2o),height,'k-');
    axis([10^0 10^4 75 110]);
    xlabel('Volume Emission Rate [cm^-3 s^-1]');
    ylabel('Altitude [km]');
    title('O2 Atmospheric Band');
end
  
    %Perturbed Emissions
    clear visOH; clear visOI;
    k6N2=5.7E-34.*(T./300).^(-2.62);
    k6O2=5.96E-34.*(T./300).^(-2.37);
    k1=4.7E-33.*(300./T).^2;
    k5=4.0E-12.*exp(-865./T);
    k1o3=1.4E-10.*exp(-470./T);
    k1o2a=(4.7e-33).*(300./T).^2;
    k2n2a=2.2e-15;
    k2o2a=4e-17;
    A1=0.079; A2=0.083;
    k1na=(1.1e-9).*exp(-116./T);
    alpha=0.093;
    A5=1.18;
    A6=1.35;
    CpO=2.1E2;
    CpO2=1.5E1;
   
    Mp=dox2+dnit2+dox;
    visOH8=.29.*dhyd.*dox3.*k1o3./(260+2E-11.*dox2);
    visOI=A5.*k1.*(dox.^3).*Mp./((A6+k5.*dox2).*(CpO2.*dox2+CpO.*dox));    
    visO2=(k1o2a.*A1.*(dox.^2).*(dox2+dnit2).*dox2)./((A2+k2o2a.*dox2+k2n2a.*dnit2).*(7.5.*dox2+33.*dox));
%    visNa=k1na.*alpha.*dna.*dox3;
   
    sk1=0.03.*k1;
    k2=5.0e-13;
    k3=6e-12;
    A1=1.0e-3;
    sk3=0.2*k3;
    nO2c=sk1.*(dox.^2).*Mp./(k2.*dox2+k3.*dox+A1);
    nO1s=sk3.*dox.*nO2c./(k5.*dox2+A6);
    visOIHickey=A5.*nO1s;
         
   
    %OH Emission (Makhlouf)
if Frame==0
    visOH3int0=trapz(height(minH:maxH).*100000,A31.*dOH30(:,:,minH:maxH),3);
    visOH3T0=(1./visOH3int0).*trapz(height(minH:maxH).*100000,A31.*dOH30(:,:,minH:maxH).*T(:,:,minH:maxH),3);
else

    visOH3int=trapz(height(minH:maxH).*100000,A31.*dOH3(:,:,minH:maxH),3);
    visOH3T=(1./visOH3int).*trapz(height(minH:maxH).*100000,A31.*dOH3(:,:,minH:maxH).*T(:,:,minH:maxH),3);
    visOH3Tlin=(1./visOH3int0).*trapz(height(minH:maxH).*100000,A31.*dOH30(:,:,minH:maxH).*T(:,:,minH:maxH),3);
    visOH3rel=(visOH3int-visOH3int0)./visOH3int0;
    visOH3Trel=(visOH3T-visOH3T0)./visOH3T0;
    visOH3Tlinrel=(visOH3Tlin-visOH3T0)./visOH3T0;

    %OI Emission (McDade)
    visOIinto=trapz(height(minH:maxH)'.*100000,visOIo(minH:maxH));
    visOIint=trapz(height(minH:maxH).*100000,visOI(:,:,minH:maxH),3);
    % EXPAND
    visOIrel=visOIint-visOIinto;
   
    %OI Emission (Hickey)
    visOIinthso=trapz(height(minH:maxH)'.*100000,visOIhso(minH:maxH));
    visOIHint=trapz(height(minH:maxH).*100000,visOIHickey(:,:,minH:maxH),3);
    % EXPAND
    visOIHrel=(visOIHint-visOIinthso)./visOIinthso;  

    %Store Airglow Image "Data"
if exist('Tevol')
    Tevol=cat(3,Tevol, visOH3Trel);
    OHevol=cat(3,OHevol, visOH3rel);
    OIevol=cat(3,OIevol, visOIHrel);

else
    Tevol=visOH3Trel;
    OHevol=visOH3rel;
    OIevol=visOIHrel;
end;

   %Plot airglow "images" and "slices" 
   fvers=figure(201);
   subplot(1,3,1);
   Tpert=Tp(minX:maxX,sliceY,minH:maxH);%./T0(minX:maxX,sliceY,minH:maxH);
   imagesc(x(minX:maxX)./1000,z(minH:maxH)./1000,squeeze(Tpert)');
   axis xy;
   axis image;
   axis([0 20 75 115]);
   title('T Pert. (K)');      
   colormap(gray(4096));
   caxis([-max(max(abs(Tpert))) max(max(abs(Tpert)))]);
   caxis([-100 100]);
   colorbar('horiz');
   subplot(1,3,2);
   imagesc(x(minX:maxX)./1000,z(minH:maxH)./1000,squeeze(visOIHickey(minX:maxX,sliceY,minH:maxH))');
   axis xy;
   axis image;
   axis([0 20 75 115]);
   title('OI (cm^{-3} s^{-1})');
   caxis([0 100]);
   colorbar('horiz');
   subplot(1,3,3);
   verOH3=A31.*((dOH3));
   imagesc(x(minX:maxX)./1000,z(minH:maxH)./1000,squeeze(verOH3(minX:maxX,sliceY,minH:maxH))');
   axis xy;
   axis image;
   axis([0 20 75 115]);
   title('OH (cm^{-3} s^{-1})');
   caxis([0 40000]);
   colorbar('horiz');
   n1 = Frame+1000;
   pfname = ['Vers',num2str(n1),'.png'];
   pfname(5) = '0';    
   print(fvers,pfname,'-dpng','-r300'); 
   
   figure(202);
   plot(x((1+minX):maxX),visOH3Tlinrel((1+minX):maxX,sliceY).*100,'k:',x((1+minX):maxX),visOH3Trel((1+minX):maxX,sliceY).*100,'k-',x((1+minX):maxX),visOH3rel((1+minX):maxX,sliceY).*100,'k--',x((1+minX):maxX),(visOH3Trel((1+minX):maxX,sliceY)-visOH3Tlinrel((1+minX):maxX,sliceY)).*100,'b-'); 
   legend('Linear BWT','BWT','Intensity','Difference');
   title('Intensity and Brightness-Weighted Temperature Perturbations');
   ylabel('Percent Perturbation (%)');
   xlabel('Streamwise Distance');
 
%    figure(203);
%    imagesc(x(minX:maxX),y,visOH3Trel'.*100);
%    axis xy;
%    axis image;
%    title('Temperature Perturbation (%)');      
%    colormap(jet(4096));
   caxis([-max(max(abs(visOH3Trel.*100))) max(max(abs(visOH3Trel.*100)))]);
%    caxis([-6 6]);
%    colorbar;
   
   fimge=figure(204);
   subplot(2,1,1);
   imagesc(x(minX:maxX),y,visOIHrel'.*100);
   axis xy;
   axis image;
   axis off;
   title('OI Perturbation (%)');      
   colormap(gray(1024));
%   caxis([-max(max(abs(visOH3Trel.*100))) max(max(abs(visOH3Trel.*100)))]);
%   caxis([-30 30]);
   colorbar;   
   subplot(2,1,2);
   imagesc(x(minX:maxX),y,visOH3rel'.*100);
   axis xy;
   axis image;
   axis off;
   title('OH(3,1) Perturbation (%)');
   colormap(gray(1024));
   caxis([-max(max(abs(visOH3rel.*100))) max(max(abs(visOH3rel.*100)))]);
%   caxis([-60 60]);
   colorbar;
%    figure(2);
    n1 = Frame+1000;
    pfname = ['Imge',num2str(n1),'.png'];
    pfname(5) = '0';
   print(fimge,pfname,'-dpng','-r300');   

   if (Frame>=MaxFrames)
      save 'Tevol.mat' Tevol;
      save 'OHevol.mat' OHevol;
      save 'OIevol.mat' OIevol;
   end;
   
end;

