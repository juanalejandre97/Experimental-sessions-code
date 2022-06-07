%% PRACTICA DE RADIACION GAMMA %%
close all

excel = 'gamma.xlsx';
rango1 = 'E7:E518';
rango2 = 'F7:F518';
rango3 = 'G7:G518';
rango4 = 'H7:H518';
rango5 = 'I7:I518';
rango6 = 'J7:J518';

rango7 = 'M7:M518';
rango8 = 'N7:N518'; 
rango9 = 'O7:O518';
rango10 = 'P7:P518';
rango11 = 'Q7:Q518';

co57 = xlsread(excel,rango1);
co60 = xlsread(excel,rango2);
cs137 = xlsread(excel,rango3);
na22 = xlsread(excel,rango4);
problema = xlsread(excel,rango5);
fondo = xlsread(excel,rango6);
csa = xlsread(excel,rango7)./180;
csb = xlsread(excel,rango8)./180;
csc = xlsread(excel,rango9)./180;
csd = xlsread(excel,rango10)./180;
cssolo = xlsread(excel,rango11)./180;
z = linspace(1,511,511);
canal = transpose(z);


%% VALORES MUESTRAS %%
smcs137 = 30.2;
smco60 = 5.3;
smco57 = 271.8;
smna22 = 2.6;

%% CALIBRACION EN ENERGIA %%

e1co60 = 1173.238;
e2co60 = 1332.502;
ecs137 = 661.660;
e2co57 = (136.4743+122.0614)/2;
e1na22 = 511.003;
e2na22 = 1274.542;

canale1co60 = 397;
canale2co60 = 450;
canalecs137 = 230;
canaleco57 = 47;
canale1na22 = 180;
canale2na22 = 432;

e = [e1co60,e2co60,ecs137,e2co57,e1na22,e2na22];
canales = [canale1co60,canale2co60,canalecs137,canaleco57,canale1na22,canale2na22];

% Analisis de regresion lineal%
m=1;   % grado del ajuste 
[p,S] = polyfit(canales,e,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,canales,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(canales,e);
fprintf('\n CALIBRACIÓN EN ENERGÍA')
fprintf('\nLos coeficientes de ajuste p(1)*x+p(2) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)
fprintf('\n La pendiente es %f +- %f',p(1),sigma(1))
fprintf('\n La ordenada en el origen es %f +- %f',p(2),sigma(2))

% CONVERSION %
energia = p(1).*canal+p(2);
errorenergia = @(x) x.*(sigma(1))+p(1)+sigma(2)
erroresenergia = errorenergia(canales);
 
%  energia = p(1).*canal+p(2);
%  errorenergiavector = energia*(((sigma(1)/canal)+(1/p(1)))+(sigma(2)/p(2)));
%   errorenergia = @(x) x.*(((sigma(1)./x)+(1/p(1)))+(sigma(2)/p(2)));
%   
%   erroresenergia = errorenergia(canales);
 
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

%% CALIBRACION EN RESOLUCION %%
load('co57gau.mat');
load('co60gau1.mat');
load('co60gau2.mat');
load('cs137gau.mat');
load('na22gau1.mat');
load('na22gau2.mat');
load('fondoco57');
load('fondo1co60');
load('fondo2co60');
load('fondocs137');
load('fondo1na22');
load('fondo2na22');
load('co57ef');
load('co60ef1');
load('co60ef2');
load('cs137ef');
load('na22ef1');
load('na22ef2');
agco57 = 240.2;
bgco57 = 119.6;
cgco57 = 12.27;
sigmaco57 = cgco57/sqrt(2);
ag1co60 = 540.6;
bg1co60 = 1168;
cg1co60 = 55.75;
sigma1co60 = cg1co60/sqrt(2);
ag2co60 = 428.8;
bg2co60 = 1326;
cg2co60 = 52.1;
sigma2co60 = cg2co60/sqrt(2);
agcs137 = 2306;
bgcs137 = 665.6;
cgcs137 = 33.93;
sigmacs137 = cgcs137/sqrt(2);
ag1na22 = 1792;
bg1na22 = 515;
cg1na22 = 32.05;
sigma1na22 = cg1na22/sqrt(2);
ag2na22 = 215.3;
bg2na22 = 1267;
cg2na22 = 49.72;
sigma2na22 = cg2na22/sqrt(2);

 fwhmcuadrado = (2*sqrt(2*log(2)))^2.*[sigma1co60^2,sigma2co60^2,sigmacs137^2,sigmaco57^2,sigma1na22^2,sigma2na22^2];
 
 % Analisis de regresion lineal%
m=1;   % grado del ajuste 
[p,S] = polyfit(e,fwhmcuadrado,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,e,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(e,fwhmcuadrado);
fprintf('\n \n CALIBRACIÓN EN RESOLUCIÓN')
fprintf('\nLos coeficientes de ajuste p(1)*x+p(2) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)
fprintf('\n La pendiente es %f +- %f',p(1),sigma(1))
fprintf('\n La ordenada en el origen es %f +- %f',p(2),sigma(2))


fwhmcuadradolin = p(1).*energia+p(2);
errorfwhmcuadrado = @(x) x.*(sigma(1))+p(1)+sigma(2);
erroresfwhmcuadradolin = errorfwhmcuadrado(e);

%errorfwhmcuadradolin = @(x,y) x.*(((sigma(1)./x)+(y/p(1)))+(sigma(2)/p(2)));
%error3 = errorfwhmcuadradolin(fwhmcuadrado,erroresenergia);

% REPRESENTACION %
 figure;
plot (e,fwhmcuadrado,'b.','MarkerSize',16)
hold on
plot (energia,fwhmcuadradolin,'r-','MarkerSize',16)
errorbar(e,fwhmcuadrado,erroresfwhmcuadradolin,'ko')
title('\fontsize{22} \bf \color{blue} Calibración en resolución')
axis tight;
legend( '\fontsize{20} Valores','\fontsize{20} Ajuste',  'Location', 'Northwest')
ax = gca;
ax.FontSize = 19;
xlabel('\fontsize{20} Energía (keV)');
ylabel('\fontsize{20} (FWHM)^{2} (keV)^{2}');

resolu = @(x) sqrt(p(1).*x+p(2))./x;
resolucion = resolu(e);

 % Analisis de regresion lineal%
m=2;   % grado del ajuste 
[p,S] = polyfit(e,resolucion,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,e,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(e,resolucion);
fprintf('\n \n CALIBRACIÓN DE LA RESOLUCIÓN')
fprintf('\nLos coeficientes de ajuste p(1)*x+p(2) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)
fprintf('\n La pendiente es %f +- %f',p(1),sigma(1))
fprintf('\n La ordenada en el origen es %f +- %f',p(2),sigma(2))

resolufinal = @(x) p(1).*x.^2+p(2).*x+p(3);
s = linspace(0,1400,1400);
errorresolufin = @(x)sigma(1)*x.^2+sigma(2)*x+sigma(3);
errorresolufinal = errorresolufin(e);


% REPRESENTACION %
 figure;
plot (e,resolucion,'b.','MarkerSize',16)
hold on
plot (s,resolufinal(s),'r-','MarkerSize',16)
errorbar(e,resolucion,errorresolufinal,'ko')
title('\fontsize{22} \bf \color{blue} Resolución frente a energía')
axis tight;
legend( '\fontsize{20} Valores','\fontsize{20} Ajuste',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 19;
xlabel('\fontsize{20} Energía (keV)');
ylabel('\fontsize{20} Resolución (%)');


%% CALIBRACION EN EFICIENCIA %%
rango1n = 'V7:V518';
rango2n = 'W7:W518';
rango3n = 'X7:X518';
rango4n = 'Y7:Y518';

co57n = xlsread(excel,rango1n);
co60n = xlsread(excel,rango2n);
cs137n = xlsread(excel,rango3n);
na22n = xlsread(excel,rango4n);

sumaco57 = sum(co57ef);
suma1co60 = sum(co60ef1);
suma2co60 = sum(co60ef2);
sumacs137 = sum(cs137ef);
suma1na22 = sum(na22ef1);
suma2na22 = sum(na22ef2);

% ncogau57 = (1/180).*co57gau;
 %nco60gau1 = (1/180).*co60gau1;
 %nco60gau2 = (1/180).*co60gau2;
% ncs137gau = (1/180).*cs137gau;
% nna22gau1 = (1/180).*na22gau1;
% nna22gau2 = (1/180).*na22gau2;
 %nfondo = (1/180).*fondo;
% nfondoco57 = (1/180).*fondoco57;
 %nfondo1co60 = (1/180).*fondo1co60;
 %nfondo2co60 = (1/180).*fondo2co60;
% nfondocs137 = (1/180).*fondocs137;
% nfondo1na22 = (1/180).*fondo1na22;
% nfondo2na22 = (1/180).*fondo2na22;

% sumaco57 = sum(ncogau57);
 %suma1co60 = sum(nco60gau1);
 %suma2co60 = sum(nco60gau2);
% sumacs137 = sum(ncs137gau);
% suma1na22 = sum(nna22gau1);
% suma2na22 = sum(nna22gau2);
% sumafondo = sum(nfondo);
% sumafondoco57 = sum(nfondoco57);
 %sumafondo1co60 = sum(nfondo1co60);
 %sumafondo2co60 = sum(nfondo2co60);
% sumafondocs137 = sum(nfondocs137);
% sumafondo1na22 = sum(nfondo1na22);
% sumafondo2na22 = sum(nfondo2na22);

% Actividades inciales
actinicial = 3.7E+4;
t0 = (365*5+127)*24*60*60;
lamco57 = log(2)/(271.8*24*60*60);
lamco60 = log(2)/(5.3*365*24*60*60);
lamcs137 = log(2)/(30.2*365*24*60*60);
lamna22 = log(2)/(2.6*365*24*60*60);

% Actividades
aco57 = actinicial*exp(-lamco57*t0);
aco60 = actinicial*exp(-lamco60*t0);
acs137 = actinicial*exp(-lamcs137*t0);
ana22 = actinicial*exp(-lamna22*t0);

% Intensidades de emision
ico57 = 0.8568;
i1co60 = 0.9989;
i2co60 = 0.9998;
ics137 = 0.847;
%i1na22 = 1.8;
i1na22 = 0.898;
i2na22 = 0.9993;

efico57 = (sumaco57)/(ico57*aco57);
 efi1co60 = (suma1co60)/(i1co60*aco60);
 efi2co60 = (suma2co60)/(i2co60*aco60);
 efics137 = (sumacs137)/(ics137*acs137);
 efi1na22 = (suma1na22)/(i1na22*ana22);
 efi2na22 = (suma2na22)/(i2na22*ana22);

% efico57 = (sumaco57-sumafondoco57)/(ico57*aco57);
 %efi1co60 = (suma1co60-sumafondo1co60)/(i1co60*aco60);
 %efi2co60 = (suma2co60-sumafondo2co60)/(i2co60*aco60);
% efics137 = (sumacs137-sumafondocs137)/(ics137*acs137);
% efi1na22 = (suma1na22-sumafondo1na22)/(i1na22*ana22);
% efi2na22 = (suma2na22-sumafondo2na22)/(i2na22*ana22);
%% AJUSTE DE LA CALIBRACION EN EFICIENCIA Y REPRESENTACION %% 
sumas = [sumaco57,suma1co60,suma2co60,sumacs137,suma1na22,suma2na22];
eficiencias = [efi1co60,efi2co60,efics137,efico57,efi1na22,efi2na22];

efilog = log(eficiencias);
elog = log(e);

% Analisis de regresion cuadratica %
m=2;   % grado del ajuste 
[p,S] = polyfit(elog,efilog,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,elog,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(elog,efilog);
fprintf('\n \n AJUSTE CUADRÁTICO')
fprintf('\n Los coeficientes de ajuste P(1)*X^2+p(2)*x+p(3) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)


a = -0.3907;
b = 3.616;
c = -11.36;
deltaamin = -0.9724;
deltaamax = 0.1911;
deltabmin = -3.428;
deltabmax = 10.66;
deltacmin = -32.26;
deltacmax = 9.549;

fwhmcuadradolin = p(1).*energia+p(2);
errorfwhmcuadrado = @(x) x.*(sigma(1))+p(1)+sigma(2);
erroresfwhmcuadradolin = errorfwhmcuadrado(e);


f = @(x)c+b*x+a*x.^2;
s = linspace(4.5,7.5,1500);

errorf = @(x)sigma(1)*x.^2+sigma(2)*x+sigma(3);
errorap = errorf(efilog);
% fmin = @(x)deltacmin+deltabmin*x+deltaamin*x.^2;
% fmax = @(x)deltacmax+deltabmax*x+deltaamax*x.^2;
% errormin = abs(f(elog)-fmin(elog));
% errormax = abs(f(elog)-fmax(elog));

figure
plot(elog,efilog,'k.','MarkerSize',16)
hold on
plot(s,f(s),'b-','MarkerSize',16)
errorbar(elog,efilog,errorap,'o')
title('\fontsize{22} \bf \color{blue} Calibración en eficiencia')
axis tight;
legend( '\fontsize{20} Valores','\fontsize{20} Ajuste',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 19;
xlabel('\fontsize{20} ln(E)');
ylabel('\fontsize{20} ln(\epsilon)');
%% COEFICIENTE DE ATENUACIÓN DEL Pb EN Cs %%
rhoxa = 1.13;
rhoxb = 3.4;
rhoxc = 6.8;
rhoxd = 10.77;

load('csan');
load('csbn');
load('cscn');
load('csdn');
load('cssolon');

sumacsa = sum(csan);
sumacsb = sum(csbn);
sumacsc = sum(cscn);
sumacsd = sum(csdn);
sumacssolo = sum(cssolon);

sumas = [sumacsa/sumacssolo,sumacsb/sumacssolo,sumacsc/sumacssolo,sumacsd/sumacssolo,];

y = -log(sumas);
x = [rhoxa,rhoxb,rhoxc,rhoxd];

% Analisis de regresion lineal%
m=1;   % grado del ajuste 
[p,S] = polyfit(x,y,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,x,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(x,y);
fprintf('\n \n COEFICIENTE DE ATENUACIÓN DEL Pb')
fprintf('\nLos coeficientes de ajuste p(1)*x+p(2) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)
fprintf('\n El coeficiente de atenuación del Pb es %f +- %f',p(1),sigma(1))
fprintf('\n La ordenada en el origen es %f +- %f',p(2),sigma(2))

error2 = @(x) sigma(1)*x+sigma(2);
errorap2 = error2(y);



% REPRESENTACION %
 figure;
plot (x,y,'b.','MarkerSize',16)
hold on
plot (x,yfit,'r-','MarkerSize',16)
errorbar(x,y,errorap2,'o')
title('\fontsize{22} \bf \color{blue} Coeficiente de atenuación del Pb')
axis tight;
legend( '\fontsize{20} Valores','\fontsize{20} Ajuste',  'Location', 'Southeast')
ax = gca;
ax.FontSize = 19;
ylabel('\fontsize{20} ln(N/N_{0})');
xlabel('\fontsize{20} \rho·x (g/cm^{2})');


 %% REPRESENTACIONES DE ESPECTROS %%
% 
 %% Co57 %%
% % Representacion %
load('co57gau');
load('eco57gau');
gauco57 = @(x) agco57*exp(-((x-bgco57)/cgco57).^2)+172;
figure;
plot (energia,co57,'r-','MarkerSize',16)
hold on
plot (eco57gau,gauco57(eco57gau),'b-','MarkerSize',16)
title('\fontsize{20} \bf \color{blue} Co57')
axis tight;
legend( '\fontsize{15} Cuentas','\fontsize{15} Ajuste del fotopico',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 16;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
cco57 = [46];
errorco57 = errorenergia(cco57);
fprintf('\n \n El error de los picos del Fe57 es %f',errorco57)
%set(gca,'yscale','log')
% 
 %% Co60 %%
% % Representacion %
%load('co60gau1');
%load('co60gau2');
load('eco60gau1');
load('eco60gau2');
gau1co60 = @(x) ag1co60*exp(-((x-bg1co60)/cg1co60).^2)+46;
gau2co60 = @(x) ag2co60*exp(-((x-bg2co60)/cg2co60).^2)+25;
figure;
plot (energia,co60,'r-','MarkerSize',18)
hold on 
plot (eco60gau1,gau1co60(eco60gau1),'b-','MarkerSize',16)
plot (eco60gau2,gau2co60(eco60gau2),'g-','MarkerSize',16)
title('\fontsize{20} \bf \color{blue} Co60')
axis tight;
legend( '\fontsize{15} Cuentas','\fontsize{15} Ajuste del fotopico 1','\fontsize{15} Ajuste del fotopico 2',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
c1co60 = [398,449];
errorc1co60 = errorenergia(c1co60);
fprintf('\n \n El error de los picos del Fe57 es %f',errorc1co60)
%% Cs137 %%
% % Representacion %
%load('cs137gau');
load('ecs137gau');
gaucs137 = @(x) agcs137*exp(-((x-bgcs137)/cgcs137).^2)+34;

figure;
plot (energia,cs137,'r-','MarkerSize',18)
hold on
plot (ecs137gau,gaucs137(ecs137gau),'b-','MarkerSize',16)
title('\fontsize{20} \bf \color{blue} Cs137')
axis tight;
legend( '\fontsize{15} Cuentas','\fontsize{15} Ajuste del fotopico',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
ccs137 = [230];
errorccs137 = errorenergia(ccs137);
fprintf('\n \n El error de los picos del cs137 es %f',errorccs137)
 %% Na 22 %%
% % Representacion %
%load('na22gau1');
%load('na22gau2');
load('ena22gau1');
load('ena22gau2');

gau1na22 = @(x) ag1na22*exp(-((x-bg1na22)/cg1na22).^2)+30;
gau2na22 = @(x) ag2na22*exp(-((x-bg2na22)/cg2na22).^2)+7;
figure;
plot (energia,na22,'m-','MarkerSize',18)
hold on 
plot (ena22gau1,gau1na22(ena22gau1),'b-','MarkerSize',16)
plot (ena22gau2,gau2na22(ena22gau2),'g-','MarkerSize',16)
title('\fontsize{20} \bf \color{blue} Na22')
axis tight;
legend( '\fontsize{15} Cuentas','\fontsize{15} Ajuste del fotopico 1','\fontsize{15} Ajuste del fotopico 2','Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
cna22 = [180,432];
errorcna22 = errorenergia(cna22);
fprintf('\n \n El error de los picos del na22 es %f',errorcna22)
 %% Problema %%
 sigmapico1 = 6.179/sqrt(2);
 sigmapico2 = 10.94/sqrt(2);
 sigmapico3 = 49.2/sqrt(2);
 sigmapico4 = 1.083/sqrt(2);
% % Representacion %
figure;
plot (energia,problema,'c-','MarkerSize',18)
%plot (canal,problema,'c-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Muestra problema')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
set(gca,'yscale','log')
cpro = [13,32,108,128];
errorpro = errorenergia(cpro);
fprintf('\n \n El error de los picos de la muestra problema es %f',errorpro)
 %% Fondo %%
% % Representacion %
figure;
plot (energia,fondo,'r-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} fondo')
axis tight;
legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Nº de cuentas');
set(gca,'yscale','log')
% 
 %% csa %%
% % Representacion %
% figure;
% plot (energia,csa,'g-','MarkerSize',18)
% title('\fontsize{20} \bf \color{blue} Cs(A)')
% axis tight;
% legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
% ax = gca;
% ax.FontSize = 18;
% xlabel('\fontsize{22} Energía (keV)');
% ylabel('\fontsize{20} Nº de cuentas');
%  %% csb %%
% % % Representacion %
% figure;
% plot (energia,csb,'r-','MarkerSize',18)
% title('\fontsize{20} \bf \color{blue} Cs(B)')
% axis tight;
% legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
% ax = gca;
% ax.FontSize = 18;
% xlabel('\fontsize{22} Energía (keV)');
% ylabel('\fontsize{20} Nº de cuentas');
%  %% csc %%
% % % Representacion %
% figure;
% plot (energia,csc,'b-','MarkerSize',18)
% title('\fontsize{20} \bf \color{blue} Cs(C)')
% axis tight;
% legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
% ax = gca;
% ax.FontSize = 18;
% xlabel('\fontsize{22} Energía (keV)');
% ylabel('\fontsize{20} Nº de cuentas');
%  %% csd %%
% % % Representacion %
% figure;
% plot (energia,csd,'k-','MarkerSize',18)
% title('\fontsize{20} \bf \color{blue} Cs(D)')
% axis tight;
% legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
% ax = gca;
% ax.FontSize = 18;
% xlabel('\fontsize{22} Energía (keV)');
% ylabel('\fontsize{20} Nº de cuentas');
%  %% cssolo %%
% % % Representacion %
% figure;
% plot (energia,cssolo,'m-','MarkerSize',18)
% title('\fontsize{20} \bf \color{blue} Cs (SOLO)')
% axis tight;
% legend( '\fontsize{15} Cuentas',  'Location', 'Northeast')
% ax = gca;
% ax.FontSize = 18;
% xlabel('\fontsize{22} Energía (keV)');
% ylabel('\fontsize{20} Nº de cuentas');

%% Cs conjunta %%
figure
plot (energia,csa,'b-','MarkerSize',18)
hold on 
plot (energia,csb,'g-','MarkerSize',18)
plot (energia,csc,'r-','MarkerSize',18)
plot (energia,csd,'k-','MarkerSize',18)
plot (energia,cssolo,'m-','MarkerSize',18)
title('\fontsize{20} \bf \color{blue} Muestra de Cs^{137} con distintas láminas de Pb')
axis tight;
legend( '\fontsize{15} Cs(A)','\fontsize{15} Cs(B)','\fontsize{15} Cs(C)','\fontsize{15} Cs(D)','\fontsize{15} Cs(Solo)',  'Location', 'Northeast')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Energía (keV)');
ylabel('\fontsize{20} Tasa de cuentas (s^{-1})');
