
S = [57.50	0.86	16.31	6.73	4.43	6.78	3.45	1.95	0.22];  % low-SiO2 mean andesite
% S = [59.98	0.85	16.48	6.25	3.11	6.00	3.22	2.29	0.20];  % mean andesite
% S = [62.64	0.78	16.41	5.21	1.90	5.14	3.27	2.75	0.17];  % high-SiO2 mean andesite

Fe_melt = 0.80*[40.48   5.59	1.89	24.35	13.48	12.10	0.68	0.34	2.26] ... % mean Fe-cpx in globules Pietruszka et al., 2022
        + 0.20*[0.41	2.51	2.50	83.41	1.88	0.17	0.02	0.01	0.00];    % mean Type-II magnetite Velasco et al., 2016
    
Si_melt =      [67.94	0.25	18.05	4.32	0.06	0.80	4.80	8.99	0.20];    % mean Si-rich dacite glass Pietruszka et al., 2022


M = [Si_melt   % high-Si dacite glass in inclusions
     Fe_melt   % Fe-globules in melt inclusions
     54.25	0.04	29.57	0.66	0.02	10.96	4.19	0.45	0.00   % plg in host andesite
     53.95	0.26	1.16	17.59	24.63	1.48	0.03	0.01	0.00   % opx in host andesite
     51.94	0.63	2.38	8.64	14.15	21.26	0.37	0.01	0.00]; % cpx in host andesite

S = S./sum(S,2)*100;
M = M./sum(M,2)*100;

C = M.'\S.';
C = C./sum(C)*100;

Sfit = C.'*M/100;
res  = S-Sfit;
resnorm = norm(res,2)./norm(S,2);

fprintf(1,'\n mixing proportions:\n\n')
fprintf(1,'   Si-rich liq.  %2.2f \n',C(1))
fprintf(1,'   Fe-rich liq.  %2.2f \n',C(2))
fprintf(1,'   plag          %2.2f \n',C(3))
fprintf(1,'   opx           %2.2f \n',C(4))
fprintf(1,'   cpx           %2.2f \n',C(5))

fprintf(1,'\n melt    : crystals  %2.2f : %2.2f',C(1)+C(2),C(3)+C(4)+C(5))

fprintf(1,'\n Si-melt : Fe-melt    %2.2f :  %2.2f\n',C(1)./(C(1)+C(2)),C(2)./(C(1)+C(2)));



