%% ALPHA %%

close all
excel1 = 'files/alpha1.xlsx';
range1 = 'B2:B1025';
range2 = 'C2:C1025';
range3 = 'D2:D1025';
range4 = 'E2:E1025';
range5 = 'F2:F1025';
range6 = 'G2:G1025';
range7 = 'H2:H1025';


i1 = xlsread(excel1,range1);
i2 = xlsread(excel1,range2);
i3 = xlsread(excel1,range3);
i4 = xlsread(excel1,range4);

iamsin = xlsread(excel1,range5);
iam1l = xlsread(excel1,range6);
iam2l = xlsread(excel1,range7);

x = zeros(1,1024);
for i=1:1024
    x(i)=i;
end
s = transpose(x);

% CALIBRACIÓN %

canal =[80,178,277,376,475,575,674,774,873,972];
energia = [1,1.99,2.99,3.98,4.97,5.97,6.97,7.97,8.97,9.97];

% Analisis de regresion lineal%
m=1;   % grado del ajuste 
[p,S] = polyfit(canal,energia,m);   % p coeficientes del polinomio de ajuste (en este caso una recta)
[yfit,~] = polyval(p,canal,S);
% Calculo de la desviacion standard de los coeficientes
iR = inv(S.R);
Cov = iR*iR';
sigma = sqrt( diag(Cov)*S.normr^2/S.df) ;
r = corrcoef(canal,energia);
fprintf('\nLos coeficientes de ajuste p(1)*x+p(2) son \n')
for i=1:m+1
    fprintf('p(%d) = %5.3f +- %5.3f\n', i, p(i), sigma(i));
end
fprintf('\n El coeficiente de correlacion para el ajuste es r^2 = %f',(r(2,1))^2)
fprintf('\n La pendiente es %f +- %f',p(1),sigma(1))


% Representacion %
figure;
plot (canal,energia,'r.','MarkerSize',16)
hold on 
plot(canal,yfit,'b-','MarkerSize',16);

errory = 0.01*ones(size(energia));
errorbar(canal,energia,errory,'ko')
errorx = ones(size(canal)); 
errorbar(canal,energia,errorx,'horizontal','ko')
title('\fontsize{20} \bf \color{blue} Recta de calibración')
axis tight;
legend( '\fontsize{15} Datos experimentales','\fontsize{15} Ajuste',  'Location', 'South')
ax = gca;
ax.FontSize = 18;
xlabel('\fontsize{22} Canales');
ylabel('\fontsize{20} Energía (MeV) ');

% PICO POLONIO %
came= 521;
eame = 5.4856;

b =eame-(p(1)*came);
fprintf('\n La ordenada en el origen es b = %f', b)

% recta de calibración %
rcal = @(x) p(1)*x+b;
sener = rcal(s); 

% REPRESENTACION EN ENERGIAS %
figure 
subplot(2,2,1)
plot(sener,i1,'r-')
subplot(2,2,2)
plot(sener,i2,'b-')
subplot(2,2,3)
plot(sener,i3,'k-')
subplot(2,2,4)
plot(sener,i4,'m-')

figure
plot(sener,i1,'r-')
title('\fontsize{20} \bf \color{blue} Espectro de radiación del ^{241}Am')
axis tight;
legend( '\fontsize{15} Datos experimentales','\fontsize{15} Ajuste',  'Location', 'Northwest')
ax = gca;
ax.FontSize = 14;
xlabel('\fontsize{16} Energía (MeV)');
ylabel('\fontsize{16} N ');

figure
plot(sener,i3,'k-')
title('\fontsize{20} \bf \color{blue} Espectro de radiación del ^{210}Po')
axis tight;
legend( '\fontsize{15} Datos experimentales','\fontsize{15} Ajuste',  'Location', 'Northwest')
ax = gca;
ax.FontSize = 14;
xlabel('\fontsize{16} Energía (MeV)');
ylabel('\fontsize{16} N ');

figure
plot(sener,i2,'b-')
title('\fontsize{20} \bf \color{blue} Espectro de radiación del ^{241}Am')
axis tight;
legend( '\fontsize{15} Datos experimentales','\fontsize{15} Ajuste',  'Location', 'Northwest')
ax = gca;
ax.FontSize = 14;
xlabel('\fontsize{16} Energía (MeV)');
ylabel('\fontsize{16} N ');

figure
plot(sener,i4,'m-','LineWidth',2)
title('\fontsize{20} \bf \color{blue} Espectro de radiación de la muestra problema')
axis tight;
legend( '\fontsize{15} Datos experimentales','\fontsize{15} Ajuste',  'Location', 'Northwest')
ax = gca;
ax.FontSize = 14;
xlabel('\fontsize{16} Energía (MeV)');
ylabel('\fontsize{16} N ');



figure
plot(sener,i2,'b-')

figure
plot(sener,i4,'m-')

figure
subplot(1,3,1)
plot(sener,iamsin,'r-')
subplot(1,3,2)
plot(sener,iam1l,'b-')
subplot(1,3,3)
plot(sener,iam2l,'k-')

figure
plot(sener,iamsin,'r-')
hold on
plot(sener,iam1l,'b-')
plot(sener,iam2l,'k-')

figure
plot(sener,iamsin,'r-')
hold on
plot(sener,iam1l,'b-')
plot(sener,iam2l,'k-')
title('\fontsize{20} \bf \color{blue} Espectro de radiación del ^{241}Am')
axis tight;
legend( '\fontsize{15} 0 láminas de Au','\fontsize{15} 1 láminas de Au','\fontsize{15} 2 láminas de Au',  'Location', 'Northwest')
ax = gca;
ax.FontSize = 14;
xlabel('\fontsize{16} Energía (MeV)');
ylabel('\fontsize{16} N ');

%% CALCULO DE ACTIVIDADES %%
r = sqrt(50E-6/pi);
d = 10E-3;

name = 199/180; 
npol = 172/180;

theta = atan(r/d);
g = (sin(theta/2))^2;

aname = name/(g*0.845);
apol = npol/g;

fprintf('\n Las actividades del americio y el polonio son respectivamente %f y %f Bq',aname,apol)
    