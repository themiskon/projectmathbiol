%% seasonal 
p.L=5; %depth meter 
p.xgrid=linspace(0,p.L,21); %  
p.tRange=1:1:365*5; %  
p.rmax=0.7; %  
p.H=30; % half saturation constant 
p.Hf=1.5; % 
p.tf=90; % 
p.K=5; % 
p.di=1/30; % 
p.dmax=3/10; % 
p.df=5/10; % 
p.tfp=250; % 
p.r=0.07; %phytoplankton metabolic cost /day (or mortality.. )
p.dp=0.1; %bckground mortality
%zooplankton parameters
p.epsz=0.7; %assimilation efficiency zooplankton
p.Ks=1; %half saturation for zooplankton μmolN/m^3
p.i0=1; %max consumption zooplankton μmolN/day
p.bz=0.004; %zooplankton clears b=ar^3, a=4*10^6, r=zoo radius (m).
p.dz=0.04; %/day
p.j=0.5; %uptake half saturation
lat=50; %latitude of lake %the model works up to 55
%fish parameters
p.epsf=0.7/10000; %fish assimilation efficiency
p.bf=0.01;% m^3/day
p.Cfmax=0.005*10000; %grams/day (μΜΝ/day)
p.preffp=0.9;%preference to phytoplankton
p.preffz=0.1;%preference to zooplankton
p.Mc= 0.8/365; %Metabolic cost /day

% p.Zoomax=2; %max zooplankton concentration

V=@(t,x) (exp(10*log(19/20)*x)-9/10).*(-0.3/2*cos(2*pi*t/365)+0.85).*(x>=0&x<=0.1)+(2/27*x-101/90).*(-0.2/2*cos(2*pi*t/365)+0.9).*(x>0.1&x<=3)+-0.9*(-0.1/2*cos(2*pi*t/365)+0.95).*(x>3); %m/day
Vf=@(t,x) 10000*-heaviside(x-1.5);
D=@(t,x) 10*exp(-1.5*x)*(0.4/2*cos(2*pi*t/365)+0.8);%m2/day
Df=@(t,x) 10*(1-heaviside(x-1.5));
f=@(t,x,P,Z,F,N,g) ((g.*N./(p.j+N))-p.r).*P-p.i0.*p.bz.*P.*Z./(p.Ks+p.bz.*P)-p.Cfmax.*F.*(p.bf.*P.*p.preffp)/(p.Cfmax+p.bf.*(p.preffp.*P+p.preffz.*Z))-p.dp.*P;
fz=@(t,x,Z,P,F) p.epsz.*p.i0.*p.bz.*P.*Z./(p.Ks+p.bz.*P)-p.dz.*Z-p.Cfmax.*F.*(p.bf.*p.preffz.*Z)/(p.Cfmax+p.bf.*(p.preffp.*P+p.preffz.*Z));
ff=@(t,x,Z,F,P) p.epsf.*p.Cfmax.*(p.bf.*P.*p.preffp+p.bf.*p.preffz.*Z)/(p.Cfmax+p.bf.*(p.preffp.*P+p.preffz.*Z)).*F-p.Mc.*F;
fn=@(t,x,P,N,g) ((p.dp)-(g.*N./(p.j+N))).*P;

C0=zeros([4*length(p.xgrid)],1);
C0(1:7,1)=linspace(0,7,7);%rand(size(xgrid));%C(t=0,xgrid) 
  C0(length(p.xgrid)+1:length(p.xgrid)+5,1)=linspace(0,5,5);
   C0(2*length(p.xgrid)+1:3*length(p.xgrid),1)=repmat(0,length(p.xgrid),1);
   C0(end-5:end,1)=15*ones(6,1);
   sensi1=floor(linspace(1,60,12));
   sensi1=[1 2 8 51];
   for sensi=1:1
[t,C] = solvePDE(fn,f,fz,ff, V,Vf, D,Df, p.tRange, p.xgrid, C0,lat,sensi); 
yy.(sprintf('case%d',sensi))=C;
   end
   save initial.mat
  %%
figure 
surface(t(1:end),-p.xgrid,(C(1:end,1:length(p.xgrid)))');
shading interp
hold on
% surf(p.tRange,p.xgrid,log10(C(:,length(p.xgrid)+1:2*length(p.xgrid)))','b');
% alpha 0.5
shading interp
title('seasonal change phytoplankton') 
xlabel('time (day)'); 
ylabel('Depth (m)'); 
zlabel('concentration (\muM m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
  %%
figure 
surface(t(1:end),-p.xgrid,C(1:end,length(p.xgrid)+1:2*length(p.xgrid))'); 
title('seasonal change zooplankton') 
xlabel('time (day)'); 
ylabel('Depth (m)'); 
zlabel('concentration (\muM m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
 
%%
figure 
surface(t(1:end),-p.xgrid,C(1:end,3*length(p.xgrid)+1:end)');
shading interp
hold on
% surf(p.tRange,p.xgrid,log10(C(:,length(p.xgrid)+1:2*length(p.xgrid)))','b');
% alpha 0.5
shading interp
title('seasonal change nutrients') 
xlabel('time (day)'); 
ylabel('Depth (m)'); 
zlabel('concentration (\muM m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
%%
figure 
surface(t(1:end),-p.xgrid,C(1:end,2*length(p.xgrid)+1:3*length(p.xgrid))'); 
title('seasonal change fish') 
xlabel('time (day)'); 
ylabel('Depth (m)'); 
zlabel('concentration (\gm^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
%% final graph 
tiledlayout(2,2)
nexttile 
surface(t(1:2*365),-p.xgrid,C(end-2*365+1:end,3*length(p.xgrid)+1:end)');
shading interp
hold on
% surf(p.tRange,p.xgrid,log10(C(:,length(p.xgrid)+1:2*length(p.xgrid)))','b');
% alpha 0.5
shading interp
title('seasonal change nutrients')  
ylabel('Depth (m)'); 
zlabel('concentration (muM m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
axis tight

nexttile
surface(t(1:2*365),-p.xgrid,(C(end-2*365+1:end,1:length(p.xgrid)))');
shading interp
hold on
% surf(p.tRange,p.xgrid,log10(C(:,length(p.xgrid)+1:2*length(p.xgrid)))','b');
% alpha 0.5
shading interp
title('seasonal change phytoplankton') 
zlabel('concentration (muM m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
axis tight

nexttile
surface(t(1:2*365),-p.xgrid,C(end-2*365+1:end,length(p.xgrid)+1:2*length(p.xgrid))'); 
title('seasonal change zooplankton') 
xlabel('time (day)'); 
ylabel('Depth (m)'); 
zlabel('concentration (muM m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
axis tight

nexttile
surface(t(1:2*365),-p.xgrid,C(end-2*365+1:end,2*length(p.xgrid)+1:3*length(p.xgrid))');
shading interp
hold on
shading interp
title('seasonal change fish') 
xlabel('time (day)'); 
zlabel('concentration (g m^-^3)'); 
shading interp 
colorbar 
set(gca,FontSize=20) 
axis tight