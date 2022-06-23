% routine to set three-phase coefficient calibration for system of 
% silicate crystals + Si-rich melt + Fe-rich melt
clear variables; close all;

% set directory paths
addpath('./ternplot');
workdir = pwd;
figsdir = pwd;

% set printing options
figname   = 'ThreePhaseCalibration';
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
axh = 8;
axw = 8;
ahs = 1.5;
avs = 1.5;
axb = 1.8;
axt = 1.2;
axl = 1.4;
axr = 1.2;
cbh = axh; cbw = 0.2;
fh = axb + 3*axh + 2*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;

% initialize figure and axes
f = figure;
set(f,UN{[1,3]},'Position',[1 10 fw fh]);
set(f,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(f,'Color','w','InvertHardcopy','off');
set(f,'Resize','off');
ax(1) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+2*axh+2*avs axw axh]);
ax(2) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+2*axh+2*avs axw axh]);
ax(3) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+2*axh+2*avs axw axh]);
ax(4) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(5) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(6) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(7) = axes(UN{[1,3]},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(8) = axes(UN{[1,3]},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(9) = axes(UN{[1,3]},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

% prepare phase proportion arrays
% phi1 = crystals; phi2 = Si-rich melt; phi3 = Fe-rich melt
np    =  200;
phi1  =  linspace(0,1,np);
phi2  =  linspace(0,1,np);
[phi2,phi1]  =  meshgrid(phi2,phi1);
phi2  =  phi2(:);  phi1  =  phi1(:);  phi3  =  1-phi2-phi1;

phi               =  [phi1,phi2,phi3];
phi(phi(:,3)<0,:) = nan;

% set pure-phase properties
eta0 = [1e+16,1e+5,1];       % pure-phase viscosities
d0   = [1e-3,1e-3,1e-3];     % characteristic size of local-scale phase constituents

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

% calculate phase permission weight functions
for i = 1:3
    Sphi(:,:,i)  =  (phi./B(i,:)).^(1./C(i,:)) ./ sum((phi./B(i,:)).^(1./C(i,:)),2);
end
for i = 1:3
    Xphi(:,:,i)  =  sum(A(i,:).*Sphi(:,:,i),2).*phi + (1-sum(A(i,:).*Sphi(:,:,i),2)).*Sphi(:,:,i);
end

% prepare for plotting ternary plots
m = 5;
grids = linspace(0,1,m+1);
grids = grids(1:end-1);
labels = num2str(grids(2:end)');
[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);
n = m-1;
[X,Y]  =  terncoords(phi(:,1),phi(:,2),phi(:,3));
tri    =  simpletri(np);

% plot top left panel
axes(ax(1));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 2;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,1,1));
caxis([0,1]); shading interp; view([0,90]);
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(1),FireAndIce);
c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[0:0.2:1]);
cpos = get(c,'position'); cpos(1) = axl; cpos(2) = axb+2*axh+1.6*avs; cpos(3) = axw; cpos(4) = 0.4;
set(c,'position',cpos,TL{:});
text(-0.12,1.00,'\textbf{(a)}~$X_\phi^{1 1}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot top middle panel
axes(ax(2));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 2;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,2,1));
caxis([0,1]); shading interp; view([0,90]);
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(2),FireAndIce);
text(-0.12,1.00,'\textbf{(b)}~$X_\phi^{1 2}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot top right panel
axes(ax(3));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 2;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,3,1));
caxis([0,1]); shading interp; view([0,90]);
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(3),FireAndIce);
text(-0.12,1.00,'\textbf{(c)}~$X_\phi^{1 3}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot top right panel
axes(ax(4));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,1,2));
caxis([0,1]); shading interp; view([0,90]);
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(4),FireAndIce);
text(-0.12,1.00,'\textbf{(d)}~$X_\phi^{2 1}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot bottom left panel
axes(ax(5));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,2,2));
caxis([0,1]); shading interp; view([0,90]); 
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]});
colormap(ax(5),FireAndIce);
c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[0:0.2:1]);
cpos = get(c,'position'); cpos(1) = axl+axw+ahs; cpos(2) = axb+axh+0.60*avs; cpos(3) = axw; cpos(4) = 0.4;
set(c,'position',cpos,TL{:});
text(-0.12,1.00,'\textbf{(e)}~$X_\phi^{2 2}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot bottom right panel
axes(ax(6));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,3,2));
caxis([0,1]); shading interp; view([0,90]); 
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(6),FireAndIce);

text(-0.12,1.00,'\textbf{(f)}~$X_\phi^{2 3}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot top right panel
axes(ax(7));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,1,3));
caxis([0,1]); shading interp; view([0,90]);
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(7),FireAndIce);
text(-0.12,1.00,'\textbf{(g)}~$X_\phi^{3 1}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot bottom left panel
axes(ax(8));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,2,3));
caxis([0,1]); shading interp; view([0,90]); 
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,3]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(8),FireAndIce);
text(-0.12,1.00,'\textbf{(h)}~$X_\phi^{3 2}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% plot bottom right panel
axes(ax(9));

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 20;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(4,1),'color',[0,0,0],LW{[1,4]},'handlevisibility','off'); 
axis equal tight;

% plot grid on ternary axes
for i = 1:n
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(n-i+2)],[y1(i+1),y2(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(n-i+2)],[y2(i+1),y3(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(n-i+2)],[y3(i+1),y1(n-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

trisurf(tri,X,Y,Xphi(:,3,3));
caxis([0,1]); shading interp; view([0,90]); 
text(x3(2:end)+0.02,y3(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,2]},VA{[1,2]});
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,4]});
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,FS{[1,3]},TX{:},HA{[1,4]},VA{[1,2]});
text( 1.02, 0.02,'$\phi^1$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
text(-0.02,-0.02,'$\phi^3$'   ,TX{:},FS{[1,4]},UN{[1,2]},HA{[1,4]},VA{[1,3]});
text(0.525,0.975,'$\phi^2$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,2]}); 
colormap(ax(9),FireAndIce);
c = colorbar('location','southoutside',UN{[1,3]},FS{[1,3]},TL{:},'Ticks',[0:0.2:1]);
cpos = get(c,'position'); cpos(1) = axl+2*axw+2*ahs; cpos(2) = axb-0.40*avs; cpos(3) = axw; cpos(4) = 0.4;
set(c,'position',cpos,TL{:});
text(-0.12,1.00,'\textbf{(i)}~$X_\phi^{3 3}$',TX{:},FS{[1,4]},UN{[1,2]},HA{[1,2]},VA{[1,3]});
hold off; drawnow;

% print figure
if printfig
    cd(figsdir);
    print(f,format,resl,rend,figname,'-loose');
    cd(workdir);
end