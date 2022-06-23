% calculate mixture of analysed El Laco mineral and melt inclusion
% compositions that reconstruct averaged host andesite bulk composition
clear all; close all; clc;

% target composition (solution)
%    SiO2   TiO2    Al2O3   FeOt    MgO     CaO     Na2O    K2O     P2O5
S = [57.50	0.86	16.31	6.73	4.43	6.78	3.45	1.95	0.22];  % low-SiO2 mean andesite (5 lowest SiO2 andesites, Tornos et al., 2017)
% S = [59.98	0.85	16.48	6.25	3.11	6.00	3.22	2.29    0.20];  % mean andesite (all andesites, Tornos et al., 2017)
% S = [62.64	0.78	16.41	5.21	1.90	5.14	3.27	2.75    0.17];  % high-SiO2 mean andesite (5 highest SiO2 andesites, Tornos et al., 2017)

%               SiO2   TiO2    Al2O3   FeOt    MgO     CaO     Na2O    K2O     P2O5
Fe_melt = 0.80*[40.48   5.59	1.89	24.35	13.48	12.10	0.68	0.34	2.26] ... % mean Fe-cpx in globules, Pietruszka et al., in review Nature Comms
        + 0.20*[0.41	2.51	2.50	83.41	1.88	0.17	0.02	0.01	0.00];    % mean Type-II magnetite, Velasco et al., 2016
    
%               SiO2   TiO2    Al2O3   FeOt    MgO     CaO     Na2O    K2O     P2O5
Si_melt =      [67.94	0.25	18.05	4.32	0.06	0.80	4.80	8.99	0.20];    % mean high-Si dacite glass, Pietruszka et al., in review Nature Comms

% mixing end-members
M = [Si_melt   % high-Si dacite glass in inclusions
     Fe_melt   % Fe-globules in melt inclusions
     54.25	0.04	29.57	0.66	0.02	10.96	4.19	0.45	0.00   % mean plg in host andesite, Velasco et al., 2016
     53.95	0.26	1.16	17.59	24.63	1.48	0.03	0.01	0.00   % mean opx in host andesite, Velasco et al., 2016
     51.94	0.63	2.38	8.64	14.15	21.26	0.37	0.01	0.00]; % mean cpx in host andesite, Velasco et al., 2016

% normalise to unit sum
S = S./sum(S,2)*100; 
M = M./sum(M,2)*100;

% least-squares fit to obtain abundances of end-members that reconstruct target composition
C = M.'\S.';
C = C./sum(C)*100;

% fitted reconstruction of target composition and residual
Sfit = C.'*M/100;
res  = S-Sfit;
resnorm = norm(res,2)./norm(S,2);

% print results
fprintf(1,'\n--target composition (major oxides)\n\n')
fprintf(1,'   SiO2   TiO2   Al2O3   FeOt    MgO     CaO     Na2O    K2O     P2O5\n')
fprintf(1,'   %2.2f  %2.2f   %2.2f   %2.2f    %2.2f    %2.2f    %2.2f   %2.2f     %2.2f\n',S(1),S(2),S(3),S(4),S(5),S(6),S(7),S(8),S(9));

fprintf(1,'\n--fitted composition (major oxides)\n\n')
fprintf(1,'   SiO2   TiO2   Al2O3   FeOt    MgO     CaO     Na2O    K2O     P2O5\n')
fprintf(1,'   %2.2f  %2.2f   %2.2f   %2.2f    %2.2f    %2.2f    %2.2f   %2.2f     %2.2f\n',Sfit(1),Sfit(2),Sfit(3),Sfit(4),Sfit(5),Sfit(6),Sfit(7),Sfit(8),Sfit(9));

fprintf(1,'\n--fitting residuals (major oxides)\n\n')
fprintf(1,'   SiO2   TiO2   Al2O3   FeOt    MgO     CaO     Na2O    K2O     P2O5\n')
fprintf(1,'   %2.2f  %2.2f   %2.2f    %2.2f   %2.2f   %2.2f    %2.2f   %2.2f    %2.2f\n',res(1),res(2),res(3),res(4),res(5),res(6),res(7),res(8),res(9));

fprintf(1,'\n mixing proportions:\n\n')
fprintf(1,'   Si-rich liq.  %2.2f \n',C(1))
fprintf(1,'   Fe-rich liq.  %2.2f \n',C(2))
fprintf(1,'   plag          %2.2f \n',C(3))
fprintf(1,'   opx           %2.2f \n',C(4))
fprintf(1,'   cpx           %2.2f \n',C(5))

fprintf(1,'\n melt    : crystals   %2.2f : %2.2f',C(1)+C(2),C(3)+C(4)+C(5))

fprintf(1,'\n Si-melt : Fe-melt     %2.2f :  %2.2f\n',C(1)./(C(1)+C(2)),C(2)./(C(1)+C(2)));



