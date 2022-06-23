% routine to fit T,C-dependent viscosity model to viscometry measurements
clear; close all;


% prep plotting parameters
MS = {'MarkerSize',10};
LW = {'LineWidth',1.5};

pur   = [0.6500  0.1500  0.5500];
grn   = [0.1500  0.7000  0.0500];
amb   = [0.8000  0.5000  0.2000];
lav   = [0.4000  0.1500  0.6000];
wht   = [1.0000  1.0000  1.0000];

c     =  0.3;
lpur  = [0.6500  0.1500  0.5500].*c + (1-c).*wht;
lgrn  = [0.1500  0.7000  0.0500].*c + (1-c).*wht;
lamb  = [0.8000  0.5000  0.2000].*c + (1-c).*wht;
llav  = [0.4000  0.1500  0.6000].*c + (1-c).*wht;


R   = 8.3145;  % universal gas constant
Tax = linspace(1000,1600,1e3).';  % T-axis for plotting


% read in analysed experimental compositions
%         SiO2       TiO2     Al2O3      FeO(t)    MgO       CaO       Na2O      K2O       P2O5
OX   =  [ 70.4500    0.3000   14.5600    2.6400    0.1700    0.3200    3.3500    7.7800         0 ;  % rhyolite (synthetic)
          64.2300    0.8700   17.2500    3.1500    1.0000    5.6100    4.5200    3.1400    0.1300 ;  % andesite (El Laco sample)
          39.5800    5.4600    1.6800   21.5700   12.8800   12.6100    0.6400    0.1800    2.0100 ;  % Fe-cpx   (synthetic)
           0.0700    0.3850    0.3000   89.3150    0.7650    0.0200         0         0         0 ]; % mt ore   (avg. 2 El Laco samples)

       
% combine major oxides into three end-members (EMs)
ind1 = [2 4 5 6 9];  % TiO2, FeOtot, MgO, CaO, P2O5 (Fe-rich EM)
ind2 = [3 7 8];      % Al2O3, Na2O, K2O (Al-Na-K-rich EM)
ind3 = [1];          % SiO2 (Silica EM)

EM  =  [sum(OX(:,ind1),2) sum(OX(:,ind2),2) sum(OX(:,ind3),2)];
EM  =  EM./sum(EM,2); % normalise to unit sum


% read in measured data
Tmp_ore = 1627;
Eta_ore = 0.18;

Tmp_cpx = [1483 1459 1435 1411 1387 1363 1339 1315];
Eta_cpx = [0.62 0.65 0.69 0.72 0.77 0.82 0.89 0.97];

Tmp_and = [ 1483  1459  1435  1411  1387  1363  1339   1315   1291   1267   1244   1220   1196    1172    1148    1124    1100];
Eta_and = [12.00 15.61 20.68 27.57 37.23 51.01 71.14 100.39 142.86 207.66 307.46 464.98 718.76 1137.84 1860.86 3127.40 5379.65];

Tmp_rhy = [ 1627  1603  1579  1555   1531   1507   1483   1459   1435   1411    1387    1363    1339    1315    1291];
Eta_rhy = [40.51 54.24 73.15 99.44 136.20 188.78 266.62 377.73 546.40 801.64 1181.70 1763.62 2677.91 4117.14 6373.63];


% plot measurements against T, EMs
figure(1); clf;
subplot(2,2,1)
semilogy(1e4./(Tmp_rhy+273.15),Eta_rhy,'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(1e4./(Tmp_and+273.15),Eta_and,'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav); 
semilogy(1e4./(Tmp_cpx+273.15),Eta_cpx,'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(1e4./(Tmp_ore+273.15),Eta_ore,'>',MS{:},LW{:},'Color',amb,'MarkerFaceColor',lamb);
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Inverse Temperature [10$^4$/K]','Interpreter','latex','FontSize',13);
ylabel('Viscosity [Pas]','Interpreter','latex','FontSize',13);
subplot(2,2,2)
semilogy(EM(1,1),Eta_rhy,'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,1),Eta_and,'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav); 
semilogy(EM(3,1),Eta_cpx,'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(4,1),Eta_ore,'>',MS{:},LW{:},'Color',amb,'MarkerFaceColor',lamb);
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Fe-rich [wt]','Interpreter','latex','FontSize',13);
subplot(2,2,3)
semilogy(EM(1,2),Eta_rhy,'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,2),Eta_and,'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav); 
semilogy(EM(3,2),Eta_cpx,'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(4,2),Eta_ore,'>',MS{:},LW{:},'Color',amb,'MarkerFaceColor',lamb);
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Al-Na-K-rich [wt]','Interpreter','latex','FontSize',13);
ylabel('Viscosity [Pas]','Interpreter','latex','FontSize',13);
subplot(2,2,4)
semilogy(EM(1,3),Eta_rhy,'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,3),Eta_and,'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav); 
semilogy(EM(3,3),Eta_cpx,'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(4,3),Eta_ore,'>',MS{:},LW{:},'Color',amb,'MarkerFaceColor',lamb);
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Silica [wt]','Interpreter','latex','FontSize',13);


% fit T-dependence of viscosity eta to obtain prefactors A0, activation energies Ea
C = [ones(length(Tmp_rhy),1),1./R./(Tmp_rhy.'+273.15)];
Eta_ft = C\log(Eta_rhy).';   A0(1) = exp(Eta_ft(1)); Ea(1) = Eta_ft(2);
C = [ones(length(Tmp_and),1),1./R./(Tmp_and.'+273.15)];
Eta_ft = C\log(Eta_and).';   A0(2) = exp(Eta_ft(1)); Ea(2) = Eta_ft(2);
C = [ones(length(Tmp_cpx),1),1./R./(Tmp_cpx.'+273.15)];
Eta_ft = C\log(Eta_cpx).';   A0(3) = exp(Eta_ft(1)); Ea(3) = Eta_ft(2);


% fit C-dependence of Ea, A0
C = [ones(3,1),EM(1:3,1),EM(1:3,2)];

Ea_ft = C\log(Ea).';
A0_ft = C\log(A0).';

Ea(4) = exp(Ea_ft(1) + Ea_ft(2).*EM(4,1) + Ea_ft(3).*EM(4,2));
A0(4) = exp(A0_ft(1) + A0_ft(2).*EM(4,1) + A0_ft(3).*EM(4,2));


% plot C-fitted A0, Ea against results from T-fitting
figure(3); clf;
subplot(2,3,1)
semilogy(EM(1,1),(Ea(1)),'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,1),(Ea(2)),'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
semilogy(EM(3,1),(Ea(3)),'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(:,1),exp(Ea_ft(1) + Ea_ft(2).*EM(:,1) + Ea_ft(3).*EM(:,2)),'k+',MS{:},LW{:});
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
ylabel('Act. Energy [J/mol]','Interpreter','latex','FontSize',13);
subplot(2,3,2)
semilogy(EM(1,2),(Ea(1)),'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,2),(Ea(2)),'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
semilogy(EM(3,2),(Ea(3)),'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(:,2),exp(Ea_ft(1) + Ea_ft(2).*EM(:,1) + Ea_ft(3).*EM(:,2)),'k+',MS{:},LW{:});
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
subplot(2,3,3)
semilogy(EM(1,3),(Ea(1)),'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,3),(Ea(2)),'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
semilogy(EM(3,3),(Ea(3)),'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(:,3),exp(Ea_ft(1) + Ea_ft(2).*EM(:,1) + Ea_ft(3).*EM(:,2)),'k+',MS{:},LW{:});
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
subplot(2,3,4)
semilogy(EM(1,1),(A0(1)),'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,1),(A0(2)),'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
semilogy(EM(3,1),(A0(3)),'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(:,1),exp(A0_ft(1) + A0_ft(2).*EM(:,1) + A0_ft(3).*EM(:,2)),'k+',MS{:},LW{:});
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Fe-rich [wt]','Interpreter','latex','FontSize',13);
ylabel('Prefactor [Pas]','Interpreter','latex','FontSize',13);
subplot(2,3,5)
semilogy(EM(1,2),(A0(1)),'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,2),(A0(2)),'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
semilogy(EM(3,2),(A0(3)),'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(:,2),exp(A0_ft(1) + A0_ft(2).*EM(:,1) + A0_ft(3).*EM(:,2)),'k+',MS{:},LW{:});
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Al-Na-K-rich [wt]','Interpreter','latex','FontSize',13);
subplot(2,3,6)
semilogy(EM(1,3),(A0(1)),'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn); hold on;
semilogy(EM(2,3),(A0(2)),'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
semilogy(EM(3,3),(A0(3)),'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
semilogy(EM(:,3),exp(A0_ft(1) + A0_ft(2).*EM(:,1) + A0_ft(3).*EM(:,2)),'k+',MS{:},LW{:});
set(gca,'TickLabelInterpreter','latex','FontSize',13,LW{:});
xlabel('Silica [wt]','Interpreter','latex','FontSize',13);


% plot fitted viscosity model against measurements
f = figure(4); clf;
set(f,'Units','centimeters','Position',[1 10 20 15]);
set(f,'PaperUnits','Centimeters','PaperPosition',[0 0 20 15],'PaperSize',[20 15]);
set(f,'Color','w','InvertHardcopy','off');
set(f,'Resize','off');

semilogy([10^4/(1000+273.15),10^4/(1000+273.15)],[10^-2,10^8],'k:',LW{:}); hold on;
semilogy([10^4/(1600+273.15),10^4/(1600+273.15)],[10^-2,10^8],'k:',LW{:});
p(1) = semilogy(1e4./(Tax+273.15),(exp(A0_ft(1) + A0_ft(2).*EM(1,1) + A0_ft(3).*EM(1,2))) .* exp(exp(Ea_ft(1) + Ea_ft(2).*EM(1,1) + Ea_ft(3).*EM(1,2))./R./(Tax+273.15)),'-',LW{:},'Color',grn,'MarkerFaceColor',lgrn);
p(2) = semilogy(1e4./(Tax+273.15),(exp(A0_ft(1) + A0_ft(2).*EM(2,1) + A0_ft(3).*EM(2,2))) .* exp(exp(Ea_ft(1) + Ea_ft(2).*EM(2,1) + Ea_ft(3).*EM(2,2))./R./(Tax+273.15)),'-',LW{:},'Color',lav,'MarkerFaceColor',llav);
p(3) = semilogy(1e4./(Tax+273.15),(exp(A0_ft(1) + A0_ft(2).*EM(3,1) + A0_ft(3).*EM(3,2))) .* exp(exp(Ea_ft(1) + Ea_ft(2).*EM(3,1) + Ea_ft(3).*EM(3,2))./R./(Tax+273.15)),'-',LW{:},'Color',pur,'MarkerFaceColor',lpur);
p(4) = semilogy(1e4./(Tax+273.15),(exp(A0_ft(1) + A0_ft(2).*EM(4,1) + A0_ft(3).*EM(4,2))) .* exp(exp(Ea_ft(1) + Ea_ft(2).*EM(4,1) + Ea_ft(3).*EM(4,2))./R./(Tax+273.15)),'-',LW{:},'Color',amb,'MarkerFaceColor',lamb);
p(5) = semilogy(1e4./(Tmp_rhy+273.15),Eta_rhy,'^',MS{:},LW{:},'Color',grn,'MarkerFaceColor',lgrn);
p(6) = semilogy(1e4./(Tmp_and+273.15),Eta_and,'s',MS{:},LW{:},'Color',lav,'MarkerFaceColor',llav);
p(7) = semilogy(1e4./(Tmp_cpx+273.15),Eta_cpx,'<',MS{:},LW{:},'Color',pur,'MarkerFaceColor',lpur);
p(8) = semilogy(1e4./(Tmp_ore+273.15),Eta_ore,'>',MS{:},LW{:},'Color',amb,'MarkerFaceColor',lamb);
axis([10^4/(1700+273.15),10^4/(962.5+273.15),0.1,10^7])
set(gca,'TickLabelInterpreter','latex','FontSize',15,LW{:});
title('Viscometry with fit','Interpreter','latex','FontSize',20);
xlabel('Inverse Temperature [10$^4$/K]','Interpreter','latex','FontSize',18);
ylabel('Viscosity [Pas]','Interpreter','latex','FontSize',18);
legend([p(5),p(6),p(7),p(8)],'Laco Si-rich incl.','Laco andesite','Laco Fe-cpx incl.','Laco mt incl.','Interpreter','latex','FontSize',14,'location','northwest');


% print final figure
print(f,'fit_viscosity_model','-dpdf','-loose')