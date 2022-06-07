%% PRACTICA DE RAYOS X %%
close all

excel = 'files/rx.xlsx';
rango1 = 'Z7:Z1030';
rango2 = 'AA7:AA1030';
rango3 = 'AB7:AB1030';
rango4 = 'AC7:AC1030';
rango5 = 'AD7:AD1030';
rango6 = 'J7:J1030';
rango7 = 'K7:K1030';


co57 = xlsread(excel,rango1);
mn54 = xlsread(excel,rango2);
cd109 = xlsread(excel,rango3);
amplata = xlsread(excel,rango4);
amcobre = xlsread(excel,rango5);
amba = xlsread(excel,rango6);
amp = xlsread(excel,rango7);
z = linspace(1,1024,1024);
canal = transpose(z);

%% CALIBRACION EN ENERGIAS %%
ekafe57 = 6.930;
ekbfe57 = 7.649;
ekamn54 = 5.898;
ekbmn54 = 6.490;
ekaag109 = 22.162;
ekb1ag109 = 24.942;
ekb2ag109 = 25.454;
ekacu = 8.047;
ekbcu = 8.904;
eame = 59.5;
elaba = 4.467;
elbba = 4.828;
ekaba = 32.191;
ekb1ba = 36.376;
ekb2ba = 37.255;

canalkafe57 = 104;
canalkbfe57 = 115;
canalkamn54 = 88;
canalkbmn54 = 97;
canalkaag109 = 359;
canalkb1ag109 = 404;
canalkb2ag109 = 413;
canalkacu = 131;
canalkbcu = 145;
canalame = 964;
canallaba = 73;
canallbba = 79;
canalkaba = 521;
canalkb1ba = 588;
canalkb2ba = 603;


e = [ekafe57,ekbfe57,ekamn54,ekbmn54,ekaag109,ekb1ag109,ekb2ag109,ekacu,ekbcu,eame,elaba,elbba,ekaba,ekb1ba,ekb2ba];
canales = [canalkafe57,canalkbfe57,canalkamn54,canalkbmn54,canalkaag109,canalkb1ag109,canalkb2ag109,canalkacu,canalkbcu,canalame,canallaba,canallbba,canalkaba,canalkb1ba,canalkb2ba];

% Analisis de regresion lineal%
m=1;   % grado del ajuste 
[p,S] = polyfit(canales,e,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,canales,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(canales,e);
fprintf('\nLos coeficientes de ajuste p(1)*x+p(2) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)
fprintf('\n La pendiente es %f +- %f',p(1),sigma(1))
fprintf('\n La ordenada en el origen es %f +- %f',p(2),sigma(2))

% CONVERSION %
 energia = p(1).*canal+p(2);
 %errorenergiavector = energia*(((sigma(1)/canal)+(1/p(1)))+(sigma(2)/p(2)));
 %errorener
 errorenergia = @(x) x.*(sigma(1))+p(1)+sigma(2);
 %errorenergia = @(x) x.*(((sigma(1)./x)+(1/p(1)))+(sigma(2)/p(2)));
 erroresenergia = errorenergia(canales);
 
 % REPRESENTACION %
 figure;
plot (canal,energia,'r-','MarkerSize',16)
hold on
plot (canales,e,'b.','MarkerSize',16)
errorbar(canales,e,erroresenergia,'ko')
title('\fontsize{32} \bf \color{blue} Calibración en energía')
axis tight;
legend( '\fontsize{20} Ajuste','\fontsize{20} Valores tabulados',  'Location', 'Southeast')
ax = gca;
ax.FontSize = 25;
xlabel('\fontsize{30} Canal');
ylabel('\fontsize{30} Energía (keV)');





%% REPRESENTACION INICIAL %%

%% Co57 %%
% Representacion %
figure;
plot (energia,co57,'r-','MarkerSize',16)
title('\fontsize{20} \bf \color{blue} Co^{57} ->  Fe^{57} ')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 16;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
cfe57 = [104,115,234];
errorfe57 = errorenergia(cfe57);
fprintf('\n \n El error de los picos del Fe57 es %f',errorfe57)
%% Mn54 %%
% Representacion %
figure;
plot (energia,mn54,'b-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Mn^{54} ->  Cr^{54}')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
cmn54 = [88,97];
errorcmn54 = errorenergia(cmn54);
fprintf('\n \n El error de los picos del  Mn54 es %f',errorcmn54)
%% Cd109 %%
% Representacion %
figure;
plot (energia,cd109,'k-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Cd^{109} ->  Ag^{109}')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
ccd109 = [359,404,412];
errorccd109 = errorenergia(ccd109);
fprintf('\n \n El error de los picos del  cd109 es %f',errorccd109)
%% Americio sobre plata %%
% Representacion %
figure;
plot (energia,amplata,'m-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Americio sobre Plata')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
set(gca,'yscale','log')
cagfluo = [284,333,358,404,413,964];
errorcagfluo = errorenergia(cagfluo);
fprintf('\n \n El error de los picos del americio sobre Ag es %f',errorcagfluo)
%% Americio sobre cobre %%
% Representacion %
figure;
plot (energia,amcobre,'c-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Muestra Americio sobre Cobre')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
set(gca,'yscale','log')
ccufluo = [104,131,144,284,408,963];
errorccufluo = errorenergia(ccufluo);
fprintf('\n \n El error de los picos del americio sobre Ag es %f',errorccufluo)
%% Americio sobre bario %%
% Representacion %
figure;
plot (energia,amba,'r-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Americio sobre Bario')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
set(gca,'yscale','log')
cbafluo = [73,79,522,588,603,963];
errorcbafluo = errorenergia(cbafluo);
fprintf('\n \n El error de los picos del americio sobre bario es %f',errorcbafluo)
%% Americio sobre muestra problema %%
% Representacion %
figure;
plot (energia,amp,'b-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Americio sobre muestra problema')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
set(gca,'yscale','log')
cp = [217,242,282,357,404,520,589,963];
errorcp = errorenergia(cp);
fprintf('\n \n El error de los picos del americio sobre bario es %f',errorcp)
%% Conjunta %%
figure;
plot (energia,co57,'r-','MarkerSize',18)
hold on
plot (energia,mn54,'b-','MarkerSize',18)
plot (energia,cd109,'k-','MarkerSize',18)
plot (energia,amplata,'g-','MarkerSize',18)
plot (energia,amcobre,'m-','MarkerSize',18)

title('\fontsize{26} \bf \color{blue} Espectros RX')
axis tight;
legend( '\fontsize{15} Co57','\fontsize{15} Mn54','\fontsize{15} Cd109','\fontsize{15} Am/Ag','\fontsize{15} Am/Cu',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 22;
xlabel('\fontsize{23} Energía (keV)');
ylabel('\fontsize{23} Nº de cuentas');