%% final graph  

tiledlayout(2,2) 

nexttile  

surface(t(1:2*365),-p.xgrid,C(end-2*365+1:end,3*length(p.xgrid)+1:end)'); 

shading interp 

hold on 

% surf(p.tRange,p.xgrid,log10(C(:,length(p.xgrid)+1:2*length(p.xgrid)))','b'); 

% alpha 0.5 

shading interp 

title('Nutrient')   

ylabel('Depth (m)');  

zlabel('Concentration (\muM m^-^3)');  

shading interp  

a=colorbar;
a.Label.String='Concentration (\muM m^-^3)';  

set(gca,FontSize=20)  

axis tight 

  

nexttile 

surface(t(1:2*365),-p.xgrid,(C(end-2*365+1:end,1:length(p.xgrid)))'); 

shading interp 

hold on 

% surf(p.tRange,p.xgrid,log10(C(:,length(p.xgrid)+1:2*length(p.xgrid)))','b'); 

% alpha 0.5 

shading interp 

title('Phytoplankton')  

zlabel('Concentration (\muM m^-^3)');  

shading interp  

b=colorbar;  
b.Label.String='Concentration (\muM m^-^3)'; 
set(gca,FontSize=20)  

axis tight 

  

nexttile 

surface(t(1:2*365),-p.xgrid,C(end-2*365+1:end,length(p.xgrid)+1:2*length(p.xgrid))');  

title('Zooplankton')  

xlabel('Time (day)');  

ylabel('Depth (m)');  

zlabel('Concentration (\muM m^-^3)');  

shading interp  

c=colorbar;  
c.Label.String='Concentration (\muM m^-^3)'; 
set(gca,FontSize=20)  

axis tight 

  

nexttile 

surface(t(1:2*365),-p.xgrid,C(end-2*365+1:end,2*length(p.xgrid)+1:3*length(p.xgrid))'); 

shading interp 

hold on 

shading interp 

title('Fish')  

xlabel('Time (day)');  

zlabel('Concentration (g m^-^3)');  

shading interp  

d=colorbar;
d.Label.String='Concentration (g m^-^3)'; 
set(gca,FontSize=20)  

axis tight 

%%
%% adv diffu

ttt=0:1:365;
xxx=(linspace(0,5,366))';
figure

tiledlayout(2,2) 
%
nexttile 

mesh(ttt,-(xxx),((exp(10*log(19/20)*xxx)-9/10).*(-0.3/2*cos(2*pi*ttt/365)+0.85).*(xxx>=0&xxx<=0.1)+(2/27*xxx-101/90).*(-0.2/2*cos(2*pi*ttt/365)+0.9).*(xxx>0.1&xxx<=3)+-0.9*(-0.1/2*cos(2*pi*ttt/365)+0.95).*(xxx>3))); 
hold on 

title('Advection (P Z)')   
xlabel('Time (day)');
ylabel('Depth (m)');  
zlabel('Advection velocity (m d^-^1)');  
e=colorbar;  
e.Label.String='Advection velocity (m d^-^1)'; 
set(gca,FontSize=20)  

axis tight 

%  
nexttile 

surf(ttt,-xxx,10*exp(-1.5*xxx).*(0.4/2*cos(2*pi*ttt/365)+0.8)); 

shading interp 

hold on 

shading interp 

title('Diffusion (N P Z)')  
xlabel('Time (day)');
ylabel('Depth (m)');
zlabel('Diffusivity (m^2 d^-^1)');  

shading interp  

f=colorbar; 
f.Label.String='Diffusivity (m^2 d^-^1)'; 
set(gca,FontSize=20)  

axis tight    

%
nexttile 

plot(-xxx,10000*-heaviside(xxx-1.5),LineWidth=5);  

title('Advection (F)')  

xlabel('Depth (m)');  
ylabel('Advection (m d^-^1)');  

shading interp   

set(gca,FontSize=20)  

axis tight 

%
nexttile 

plot(-xxx,10*(1-heaviside(xxx-1.5)),LineWidth=5); 

hold on 

title('Diffusion (F)')  
xlabel('Depth (m)');  
ylabel('Diffusivity (m^2 d^-^1)');  

set(gca,FontSize=20)  

axis tight 


