
clear variables; close all;

% set directory paths
addpath('/Users/tokeller/Documents/Software/ternplot/');
workdir = pwd;
figsdir = pwd;

% set printing options
figname   = 'ScalingAnalysis';
format    = '-dpdf';
resl      = '-r250';
rend      = '-opengl';
printfig  = true;

cd(workdir);

% load custom colormap
load FireAndIce.mat
color = [0.0000  0.4470  0.7410 ; ...
         0.8500  0.3250  0.0980 ; ...
         0.9290  0.6940  0.1250 ; ...
         0.4940  0.1840  0.5560 ; ...
         0.0000  0.0000  0.0000];
     
% prepare formating options
HA = {'HorizontalAlignment','left','center','right'};
VA = {'VerticalAlignment','bottom','middle','top'};
UN = {'Units','Normalized','Centimeters'};
TX = {'Interpreter','Latex'};
TL = {'TickLabelInterpreter','Latex'};
LW = {'LineWidth',1,2,3};
FS = {'FontSize',13,16,19,24,28};
MS = {'MarkerSize',6,8,12};
LS = {'LineStyle','-','--','-.',':'};
LC = {'Color',color};

% prepare axes/borders dimensions
axh = 6;
axw = 9;
ahs = 2.5;
avs = 2.8;
axb = 1.6;
axt = 0.8;
axl = 1.4;
axr = 1.4;
cbh = axh; cbw = 0.2;
fh = axb + 1*axh + 0*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;

% initialize figure and axes
f = figure;
set(f,UN{[1,3]},'Position',[1 10 fw fh]);
set(f,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off');
set(f,'Resize','off');
ax(1) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(2) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(3) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

% compute regime transition example
np    =  600;
phi1  =  linspace(0,1,np);
phi2  =  linspace(0,1,np);
[phi2,phi1]  =  meshgrid(phi2,phi1);
phi2  =  phi2(:);  phi1  =  phi1(:);  phi3  =  1-phi2-phi1;

phi               =  [phi1,phi2,phi3];
phi(phi(:,3)<0,:) = nan;

%*****  set pure phase properties  ****************************************
g0   = 9.81;                 % gravity in vertical and horizontal direction
rho0 = [2700, 2400, 4000];   % pure-phase densities
eta0 = [1e+16,1e+5,1];       % pure-phase viscosities
d0   = [1e-3 ,1e-3,1e-3];    % characteristic size of local-scale phase constituents

kv   = eta0;           % set momentum diffusivity
kf   = d0.^2./eta0;    % set volume diffusivity
Mv   = kv./kv.';       % get momentum diffusivity contrasts
Mf   = kf./kf.';       % get volume diffusivity contrasts

% set fitting parameters
A  =  [ 0.25, 0.25, 0.25; ...
        0.25, 0.25, 0.25; ...
        0.25, 0.25, 0.25; ];  % permission slopes
    
B  =  [ 0.44, 0.18, 0.38; ...
        0.60, 0.02, 0.38; ...
        0.72, 0.25, 0.03; ];  % permission step locations
    
C  =  [ 0.30, 0.30, 0.30; ...
        0.60, 0.60, 0.12; ...
        0.60, 0.12, 0.60; ];  % permission step widths

% calculate phase permission weights functions
for i = 1:3
    Sphi(:,:,i)  =  (phi./B(i,:)).^(1./C(i,:)) ./ sum((phi./B(i,:)).^(1./C(i,:)),2);
end
for i = 1:3
    Xphi(:,:,i)  =  sum(A(i,:).*Sphi(:,:,i),2).*phi + (1-sum(A(i,:).*Sphi(:,:,i),2)).*Sphi(:,:,i);
end

% calculate phase momentum permission functions
for i = 1:3
    thtv(:,i)    =  prod(Mv(i,:).^squeeze(Xphi(:,:,i)),2);
end

% calculate phase volume permission functions
for i = 1:3
    thtf(:,i)    =  prod(Mf(i,:).^squeeze(Xphi(:,:,i)),2);
end

% get phase flux coefficients
Kvi = phi.*kv.*thtv;  % momentum
Kfi = phi.*kf.*thtf;  % volume

% get phase transfer coefficients
Cvi = phi.*(1-phi).*kv./d0.^2.*thtv;  % momentum
Cfi = phi.*(1-phi).*kf./d0.^2.*thtf;  % volume

% get coefficient-based weights
ovi = Cvi./sum(Cvi,2);
ofi = Cfi./sum(Cfi,2);

% get phase segregation coefficients
Cseg = phi.^2./Cvi;

% get phase compaction coefficients
Ccmp = phi.^2./Cfi;
    
deltals = sqrt(Cseg(:,2).*Ccmp(:,1));
deltasl = sqrt(Cseg(:,1).*Ccmp(:,2));

Drho0   = abs(rho0 - sum(phi.*rho0,2));
l0      = 1;
us0     = Cseg.*Drho0.*g0;                      % segregation speed
uc0     = 0.01.*Drho0.*9.81.*l0^2./sum(Kvi,2);  % convective speed

% get segregation-compaction length scales
delta1  =  sqrt(Cseg(:,1) .* Ccmp);
delta2  =  sqrt(Cseg(:,2) .* Ccmp);
delta3  =  sqrt(Cseg(:,3) .* Ccmp);


% prepare for plotting ternary plots
m        = 5;
grids    = linspace(0,1,m+1);
grids    = grids(1:end-1);
labels   = num2str(grids(2:end)');
[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);
n        = m-1;
[X,Y]    = terncoords(phi(:,2),phi(:,3),phi(:,1));
tri      = simpletri(np);


% plot top right panel
axes(ax(1));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0.2,0,1,0.8],[sin(10/9/pi),0,0,sin(10/9/pi)],maxz.*ones(4,1),'k-',LW{[1,4]},'handlevisibility','off'); 
plot3([0.8,0.2],[sin(10/9/pi),sin(10/9/pi)],maxz.*ones(2,1),'k:',LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
plot3([0.1,0.9],sin(10/9/pi)/2.*[1,1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.1,0.9],sin(10/9/pi)/2.*[1,1],maxz.*[1,1],'k:',LW{[1,2]});

plot3([0.2,0.1],[0,sin(10/9/pi)/2],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.2,0.1],[0,sin(10/9/pi)/2],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.4,0.2],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.4,0.2],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.6,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.6,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.8,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.8,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 

plot3([0.8,0.9],[0,sin(10/9/pi)/2],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.8,0.9],[0,sin(10/9/pi)/2],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.6,0.8],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.6,0.8],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.4,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.4,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.2,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.2,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 

% plot path
plot3([0.8,0.45,0.2],[0,sin(10/9/pi)/4,sin(10/9/pi)/1],maxz.*ones(3,1),'w'  ,LW{[1,3]}); 
plot3([0.8,0.45,0.2],[0,sin(10/9/pi)/4,sin(10/9/pi)/1],maxz.*ones(3,1),'k--',LW{[1,3]}); 

uc0(phi(:,3)>0.4,:) = nan;
trisurf(tri,X,Y,log10(uc0(:,3)));
caxis([log10(1e-3./3600/24/365.25),log10(1./3600)]); shading interp; view([0,90]); 
text(x3(2:3)+0.02,y3(2:3)-0.02,labels([4 3],:),FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:5)-0.00,y2(2:5)-0.02,flipud(labels),FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(4:5)-0.02,y1(4:5)-0.02,labels([2 1],:),FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text(x2([4 3])+0.00,(y3(3)+0.01).*[1,1],labels([2 1],:),FS{[1,3]},TX{:},HA{[1,3]},VA{[1,2]});
text( 1.01,-0.02,'$\phi^\mathrm{mSi}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]});
text(-0.02,-0.06,'$\phi^\mathrm{xtl}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,4]});
text( 0.34, 1.05,'$\phi^\mathrm{mFe}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,2]});
colormap(ax(1),(FireAndIce));
c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[log10(1e-3./3600/24/365.25),log10(1./3600/24/365.25),log10(1./3600),0],'Ticklabels',{'mm/yr';'m/yr';'m/hr';'m/s'});
cpos = get(c,'position'); cpos(1) = axl; cpos(2) = axb-0.15*avs; cpos(3) = axw; cpos(4) = 0.4;
set(c,'position',cpos,TL{:});
text(-0.10,1.50,'\textbf{(a)}~magma/mush convection speed',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;


% plot bottom left panel
axes(ax(2));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0.2,0,1,0.8],[sin(10/9/pi),0,0,sin(10/9/pi)],maxz.*ones(4,1),'k-',LW{[1,4]},'handlevisibility','off'); 
plot3([0.8,0.2],[sin(10/9/pi),sin(10/9/pi)],maxz.*ones(2,1),'k:',LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
plot3([0.1,0.9],sin(10/9/pi)/2.*[1,1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.1,0.9],sin(10/9/pi)/2.*[1,1],maxz.*[1,1],'k:',LW{[1,2]});

plot3([0.2,0.1],[0,sin(10/9/pi)/2],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.2,0.1],[0,sin(10/9/pi)/2],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.4,0.2],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.4,0.2],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.6,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.6,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.8,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.8,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 

plot3([0.8,0.9],[0,sin(10/9/pi)/2],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.8,0.9],[0,sin(10/9/pi)/2],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.6,0.8],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.6,0.8],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.4,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.4,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.2,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.2,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 

% plot path
plot3([0.8,0.45,0.2],[0,sin(10/9/pi)/4,sin(10/9/pi)/1],maxz.*ones(3,1),'w'  ,LW{[1,3]}); 
plot3([0.8,0.45,0.2],[0,sin(10/9/pi)/4,sin(10/9/pi)/1],maxz.*ones(3,1),'k--',LW{[1,3]}); 

us0(phi(:,3)>0.4,:) = nan;
trisurf(tri,X,Y,log10(us0(:,3)));
caxis([log10(1e-3./3600/24/365.25),log10(1./3600)]); shading interp; view([0,90]); 
text(x3(2:3)+0.02,y3(2:3)-0.02,labels([4 3],:),FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:5)-0.00,y2(2:5)-0.02,flipud(labels),FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(4:5)-0.02,y1(4:5)-0.02,labels([2 1],:),FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text(x2([4 3])+0.00,(y3(3)+0.01).*[1,1],labels([2 1],:),FS{[1,3]},TX{:},HA{[1,3]},VA{[1,2]});
text( 1.01,-0.02,'$\phi^\mathrm{mSi}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]});
text(-0.02,-0.06,'$\phi^\mathrm{xtl}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,4]});
text( 0.34, 1.05,'$\phi^\mathrm{mFe}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,2]});
colormap(ax(2),FireAndIce);
c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[log10(1e-3./3600/24/365.25),log10(1./3600/24/365.25),log10(1./3600),0],'Ticklabels',{'mm/yr';'m/yr';'m/hr';'m/s'});
cpos = get(c,'position'); cpos(1) = axl+axw+ahs; cpos(2) = axb-0.15*avs; cpos(3) = axw; cpos(4) = 0.4;
set(c,'position',cpos,TL{:});
text(-0.10,1.5,'\textbf{(b)}~Fe-rich melt segregation speed',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;


% plot bottom right panel
axes(ax(3));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0.2,0,1,0.8],[sin(10/9/pi),0,0,sin(10/9/pi)],maxz.*ones(4,1),'k-',LW{[1,4]},'handlevisibility','off'); 
plot3([0.8,0.2],[sin(10/9/pi),sin(10/9/pi)],maxz.*ones(2,1),'k:',LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
plot3([0.1,0.9],sin(10/9/pi)/2.*[1,1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.1,0.9],sin(10/9/pi)/2.*[1,1],maxz.*[1,1],'k:',LW{[1,2]});

plot3([0.2,0.1],[0,sin(10/9/pi)/2],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.2,0.1],[0,sin(10/9/pi)/2],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.4,0.2],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.4,0.2],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.6,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.6,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.8,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.8,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 

plot3([0.8,0.9],[0,sin(10/9/pi)/2],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.8,0.9],[0,sin(10/9/pi)/2],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.6,0.8],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.6,0.8],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.4,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.4,0.6],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 
plot3([0.2,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'w' ,LW{[1,2]}); 
plot3([0.2,0.4],[0,sin(10/9/pi)/1],maxz.*[1,1],'k:',LW{[1,2]}); 

% plot path
plot3([0.8,0.45,0.2],[0,sin(10/9/pi)/4,sin(10/9/pi)/1],maxz.*ones(3,1),'w'  ,LW{[1,3]}); 
plot3([0.8,0.45,0.2],[0,sin(10/9/pi)/4,sin(10/9/pi)/1],maxz.*ones(3,1),'k--',LW{[1,3]}); 

delta3(phi(:,3)>0.4,:) = nan;
trisurf(tri,X,Y,log10(max(delta3(:,1),delta3(:,2))));
caxis([-3,3]); shading interp; view([0,90]); 
text(x3(2:3)+0.02,y3(2:3)-0.02,labels([4 3],:),FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:5)-0.00,y2(2:5)-0.02,flipud(labels),FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(4:5)-0.02,y1(4:5)-0.02,labels([2 1],:),FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text(x2([4 3])+0.00,(y3(3)+0.01).*[1,1],labels([2 1],:),FS{[1,3]},TX{:},HA{[1,3]},VA{[1,2]});
text( 1.01,-0.02,'$\phi^\mathrm{mSi}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]});
text(-0.02,-0.06,'$\phi^\mathrm{xtl}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,4]});
text( 0.34, 1.05,'$\phi^\mathrm{mFe}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,2]});
colormap(ax(3),FireAndIce);
c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[-3:3:3],'Ticklabels',{'mm';'m';'km'});
cpos = get(c,'position'); cpos(1) = axl+2*axw+2*ahs; cpos(2) = axb-0.15*avs; cpos(3) = axw; cpos(4) = 0.4;
set(c,'position',cpos,TL{:});
text(-0.10,1.5,'\textbf{(c)}~segregation length scale',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;


% print figure
if printfig
    cd(figsdir);
    print(f,format,resl,rend,figname,'-loose');
    cd(workdir);
end
