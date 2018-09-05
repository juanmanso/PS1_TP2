
config_m;


%% Parámetros de simulación
nfft = 1024;  % tamaño de la FFT 
w = (0:nfft/2-1)/nfft*pi*2; 
w = w/pi; % Normalizo
hEA = load('SEA.mat');
fs = 44100;

% Saco la data del hEA
h = hEA.h;


% Selector de ejs
bool_ej = [1 1 1];
bool_audio = 0;

if(bool_audio)
	original1=wavread('pista_1.wav');
	original2=wavread('pista_2.wav');
	original3=wavread('pista_3.wav');
	original4=wavread('pista_4.wav');
	original5=wavread('pista_5.wav');

	aud_ea1=conv(h,original1);
	aud_ea2=conv(h,original2);
	aud_ea3=conv(h,original3);
	aud_ea4=conv(h,original4);
	aud_ea5=conv(h,original5);

	wavwrite(aud_ea1,fs,'./Pistas/pista_1_EA.wav');
	wavwrite(aud_ea2,fs,'./Pistas/pista_2_EA.wav');
	wavwrite(aud_ea3,fs,'./Pistas/pista_3_EA.wav');
	wavwrite(aud_ea4,fs,'./Pistas/pista_4_EA.wav');
	wavwrite(aud_ea5,fs,'./Pistas/pista_5_EA.wav');
end






if(bool_ej(1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ANÁLISIS DEL SEA %%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grafico la transformada

H=fft(h,nfft);
H=H(1:end/2);


plot(w,abs(H),'LineWidth',4)
xlabel('w')
ylabel('|H|')
hold on

% Límites de interés

f_bottom = 20;
f_top = 16e3;

w_bottom = f_bottom/fs*2;
w_top = f_top/fs*2;

plot([w_top, w_top],[0,3.7],'color',[0 .5 0],'linewidth',7)
%plot([w_bottom, w_bottom],[0,4],'color',[0 .5 0],'linewidth',4)


% Busco los extremos locales
[PKS, LOC, EXTRA]=findpeaks(abs(H'));
plot(LOC./nfft*2,PKS,'.r','MarkerSize',10)


[PKS2, LOC2, EXTRA]=findpeaks(abs(1./H(1:f_top/fs*nfft-20)'));
plot(LOC2./nfft*2,abs(H(LOC2)),'r.','MarkerSize',10)

% Bandas de interés
LOC=sort(LOC);
LOC2=sort(LOC2);
aux=[H(LOC2(1)) H(LOC(1)) H(LOC2(2)) H(LOC(2)) H(LOC2(end))];
pks=abs(aux);

w1 = (LOC(1)+LOC2(1))/2/(nfft/2);
w2 = (LOC2(2)+LOC(1))/2/(nfft/2);
w3 = (LOC(3)+LOC2(2))/2/(nfft/2);
w4 = (LOC(end)+LOC2(4))/2/(nfft/2);

wc = [w1 w2 w3 w4];

plot([w1 w1],[0 4],'-.','color',[0.8 0.8 0],'Linewidth',7);
plot([w2 w2],[0 4],'-.','color',[0.8 0.8 0],'Linewidth',7);
plot([w3 w3],[0 4],'-.','color',[0.8 0.8 0],'Linewidth',7);
plot([w4 w4],[0 3.7],'-.','color',[0.8 0.8 0],'Linewidth',7);

plot(0:0.01:w1,(0:0.01:w1)*0 +pks(1)-0.023,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(0:0.01:w1,(0:0.01:w1)*0 +pks(1)+0.023,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w1:0.01:w2,(w1:0.01:w2)*0 +pks(2)-0.003,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w1:0.01:w2,(w1:0.01:w2)*0 +pks(2)+0.003,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w2:0.01:w3,(w2:0.01:w3)*0 +pks(3)-0.039,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w2:0.01:w3,(w2:0.01:w3)*0 +pks(3)+0.039,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w3:0.01:w4,(w3:0.01:w4)*0 +pks(4)-0.006,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w3:0.01:w4,(w3:0.01:w4)*0 +pks(4)+0.006,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w4:0.01:w_top,(w4:0.01:w_top)*0 +pks(5)-0.013,'-.','color',[0.8 0 0.8],'Linewidth',7);
plot(w4:0.01:w_top,(w4:0.01:w_top)*0 +pks(5)+0.013,'-.','color',[0.8 0 0.8],'Linewidth',7);



lab_x = ["Frecuencia normalizada $w$"];
lab_y = ["Amplitud [\\si{\\dB}]"];
leyenda = ["Respuesta en frecuencia del sistema electroacústico"];
leyenda = [leyenda; "Frecuencia máxima perceptible por el oido humano"];
%leyenda = [leyenda; "Frecuencias de límites entre bandas"];
%leyenda = [leyenda; "Tolerancias correspondientes a cada banda"];
%leyenda = [leyenda; "Máximos y mínimos obtenidos por \\emph{findpeaks}"];
loc = 'NorthEast';
AXIS = [0 1 000 04];

set_graph('plot',[lab_x; lab_y],leyenda, loc, AXIS, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DISEÑO DEL FILTRO MULTIBANDA %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
aux=[H(LOC2(1)) H(LOC(1)) H(LOC2(2)) H(LOC(2)) H(LOC2(end))];
a = abs(1./aux);

% Teórico
%	wc = [0.05 0.12 0.203 0.456]
%	a = [2.33 0.284 3.93 0.569 1.2666];

% Copia teórica de HSEA
%%%% RECTWIN %%%%
%	M = 41;
%	wc = [0.050 0.120 0.193 0.456];
%	a = aux;
%	tipo_ventana = @rectwin;
%
%%%% HAMMING %%%%
%	M = 90;
%	wc = [0.050 0.120 0.193 0.456];
%	a = aux;
%	tipo_ventana = @hamming;

%	a = [aux(1) aux(2)-0.6 aux(3) aux(4) aux(5)];
%	a = abs(1./a);

%%%% Filtro ecualizador %%%
	M = 90;
	wc = [0.028 0.140 0.178 0.456];
	a = [3.00 0.260 04.83 0.569 1.266];
	tipo_ventana = @hamming;


[heq, M] = multibanda(a,wc,M,tipo_ventana);

HEQ = fft(heq,nfft);
HEQ = HEQ(1:end/2);

figure
plot(w,abs(H),'LineWidth',4)
xlabel(lab_x)
ylabel(lab_y)
grid minor
hold on
plot(w,abs(HEQ),'r','LineWidth',4)

%	plot([w_top, w_top],[0,3.9],'color',[0 .5 0],'linewidth',7)
%	
%	plot([wc(1) wc(1)],[-00 05],'LineWidth',4,'color',[0.8 0.8 0],'.-')
%	plot(0:0.01:w1,(0:0.01:w1)*0 +pks(1)-0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot([wc(2) wc(2)],[-00 05],'LineWidth',4,'color',[0.8 0.8 0],'.-')
%	plot([wc(3) wc(3)],[-00 05],'LineWidth',4,'color',[0.8 0.8 0],'.-')
%	plot([wc(4) wc(4)],[-00 3.9],'LineWidth',4,'color',[0.8 0.8 0],'.-')
%	
%	plot(0:0.01:w1,(0:0.01:w1)*0 +pks(1)+0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w1:0.01:w2,(w1:0.01:w2)*0 +pks(2)-0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w1:0.01:w2,(w1:0.01:w2)*0 +pks(2)+0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w2:0.01:w3,(w2:0.01:w3)*0 +pks(3)-0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w2:0.01:w3,(w2:0.01:w3)*0 +pks(3)+0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w3:0.01:w4,(w3:0.01:w4)*0 +pks(4)-0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w3:0.01:w4,(w3:0.01:w4)*0 +pks(4)+0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w4:0.01:w_top,(w4:0.01:w_top)*0 +pks(5)-0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	plot(w4:0.01:w_top,(w4:0.01:w_top)*0 +pks(5)+0.04,'-.','color',[0.8 0 0.8],'Linewidth',7);
%	
%	lab_x = ["Frecuencia normalizada $w$"];
%	lab_y = ["Amplitud [\\si{\\dB}]"];
%	leyenda = ["Respuesta en frecuencia del sistema electroacústico"];
%	leyenda = [leyenda; "Seguidor con ventana de \\emph{Hamming} con M=60"];
%	leyenda = [leyenda; "Seguidor con ventana de \\emph{Hamming} con M=70"];
%	leyenda = [leyenda; "Seguidor con ventana de \\emph{Hamming} con M=80"];
%	leyenda = [leyenda; "Seguidor con ventana de \\emph{Hamming} con M=90"];
%	leyenda = [leyenda; "Frecuencia máxima perceptible por el oido humano"];
%	leyenda = [leyenda; "Frecuencias de límites entre bandas"];
%	leyenda = [leyenda; "Tolerancias correspondientes a cada banda"];
%	loc = 'NorthEast';
%	AXIS = [0 1 000 05];
%	
%	set_graph('plot',[lab_x; lab_y],leyenda, loc, AXIS, 1);
%	

%aux1=abs(fft(multibanda([a(1) 0],wc(1),M,tipo_ventana),nfft));
%plot(w,aux1(1:end/2),'color',[0 .5 0],'LineWidth',4)
%
%aux1=abs(fft(multibanda([0 a(2) 0],[wc(2) wc(1)],M,tipo_ventana),nfft));
%plot(w,aux1(1:end/2),'color',[0 .9 0],'LineWidth',4)


%%%%%%%
%% Gráfico de la rta del sistema ecualizado en dB

figure(3)
plot(w,mag2db(abs(HEQ.*H)),'LineWidth',4)
hold on
plot([w_top w_top],[-20 08],'LineWidth',4,'g')
plot(w,(w.*0)+2,'LineWidth',4,'r')
plot(w,(w.*0)-2,'LineWidth',4,'r')

lab_y2 = ["Amplitud [\\si{\\dB}]"];
leyenda = ["Respuesta en frecuencia del sistema ecualizado por FIR FLG"];
leyenda = [leyenda; "Frecuencia máxima perceptible por el oido humano"];
leyenda = [leyenda; "Tolerancia de $\\pm \\SI{2}{\\dB}$"];
loc = 'NorthEast';
AXIS = [0 1 -10 10];

set_graph('plot',[lab_x; lab_y2],leyenda, loc, AXIS, 1);

plot([wc(1) wc(1)],[-20 20],'LineWidth',4,'color',[0.8 0.8 0])
plot([wc(2) wc(2)],[-20 20],'LineWidth',4,'color',[0.8 0.8 0])
plot([wc(3) wc(3)],[-20 20],'LineWidth',4,'color',[0.8 0.8 0])
plot([wc(4) wc(4)],[-20 08],'LineWidth',4,'color',[0.8 0.8 0])


retgrp_heq = grpdelay(heq,1,nfft);
retgrp_hsea= grpdelay(h,1,nfft);
retgrp_htot= grpdelay(conv(heq,h),1,nfft);

figure
wfull=(0:nfft-1)/nfft;
hold on
plot(wfull,retgrp_heq,'LineWidth',4)
plot(wfull,retgrp_hsea,'LineWidth',4,'r')
plot(wfull,retgrp_htot,'LineWidth',4,'color',[0 0.5 0])

lab_y3 = ["Retardo de fase"];
leyenda = ["Filtro ecualizador"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema EA"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema total"];
loc = 'North';
AXIS = [0 1 40 120];

set_graph('plot',[lab_x; lab_y3],leyenda, loc, AXIS, 1);

if(bool_audio)
	wavwrite(conv(aud_ea1,heq),fs,'./Pistas/pista_1_TOT_flg.wav');
	wavwrite(conv(aud_ea2,heq),fs,'./Pistas/pista_2_TOT_flg.wav');
	wavwrite(conv(aud_ea3,heq),fs,'./Pistas/pista_3_TOT_flg.wav');
	wavwrite(conv(aud_ea4,heq),fs,'./Pistas/pista_4_TOT_flg.wav');
	wavwrite(conv(aud_ea5,heq),fs,'./Pistas/pista_5_TOT_flg.wav');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(bool_ej(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Ejercicio 1.2 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 060;
dw = 0.025;
wc = [0.045 w2 w3 w4];
wc = [wc(1)-3e-3 w2+13e-3 w3-03e-3 w4];
w_aux1= [0 wc(1)-dw]; 
w_aux2=	[wc(1)+dw wc(2)-dw]; 
w_aux3=	[wc(2)+dw wc(3)-dw]; 
w_aux4=	[wc(3)+dw wc(4)-dw]; 
w_aux5=	[wc(4)+dw w_top-dw]; 
w_aux6=	[w_top+dw 1];
w = [w_aux1 w_aux2 w_aux3 w_aux4 w_aux5 w_aux6];
a = [2.33 0.284 3.93 0.569 1.2666 1];
a = [a(1) a(1) a(2) a(2) a(3) a(3) a(4) a(4) a(5) a(5) a(6) a(6)];
W = [35. 300 051 05 001.5 10];


heq2 = firls(M,w,a,W);
HEQ2 = fft(heq2,nfft)';
HEQ2 = HEQ2(1:end/2);

% Grafico la transformada
figure
w = (0:nfft/2-1)/nfft*pi*2; 
w = w/pi;

plot(w,abs(H),'LineWidth',4)
hold on
plot(w,abs(HEQ2),'r','LineWidth',4)
grid minor
plot([wc(1) wc(1)],[-00 05],'LineWidth',4,'y')
plot([wc(2) wc(2)],[-00 05],'LineWidth',4,'y')
plot([wc(3) wc(3)],[-00 05],'LineWidth',4,'y')
plot([wc(4) wc(4)],[-00 05],'LineWidth',4,'y')

%% Gráfico de la rta del sistema ecualizado en dB
figure
plot(w,mag2db(abs(HEQ2.*H)),'LineWidth',4)
hold on
plot([w_top w_top],[-20 09],'LineWidth',4,'g')
plot(w,(w.*0)+2,'LineWidth',4,'r')
plot(w,(w.*0)-2,'LineWidth',4,'r')

lab_x = ["Frecuencia $\\frac{w}{\\pi}$"];
lab_y = ["Amplitud [\\si{\\dB}]"];
leyenda = ["Respuesta en frecuencia del sistema ecualizado por cuadrados mínimos"];
leyenda = [leyenda; "Frecuencia máxima perceptible por el oido humano"];
loc = 'NorthEast';
AXIS = [0 1 -10 10];

set_graph('plot',[lab_x; lab_y],leyenda, loc, AXIS, 1);

plot([wc(1) wc(1)],[-20 20],'LineWidth',4,'y')
plot([wc(2) wc(2)],[-20 09],'LineWidth',4,'y')
plot([wc(3) wc(3)],[-20 09],'LineWidth',4,'y')
plot([wc(4) wc(4)],[-20 09],'LineWidth',4,'y')


retgrp_heq2 = grpdelay(heq2,1,nfft);
retgrp_htot2= grpdelay(conv(heq2,h),1,nfft);

figure
wfull=(0:nfft-1)/nfft;
hold on
plot(wfull,retgrp_heq2,'LineWidth',4)
plot(wfull,retgrp_hsea,'LineWidth',4,'r')
plot(wfull,retgrp_htot2,'LineWidth',4,'color',[0 0.5 0])

lab_y3 = ["Retardo de fase"];
leyenda = ["Filtro ecualizador"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema EA"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema total"];
loc = 'North';
AXIS = [0 1 20 100];

set_graph('plot',[lab_x; lab_y3],leyenda, loc, AXIS, 1);


return;

if(bool_audio)
	wavwrite(conv(aud_ea1,heq2),fs,'./Pistas/pista_1_TOT_ls.wav');
	wavwrite(conv(aud_ea2,heq2),fs,'./Pistas/pista_2_TOT_ls.wav');
	wavwrite(conv(aud_ea3,heq2),fs,'./Pistas/pista_3_TOT_ls.wav');
	wavwrite(conv(aud_ea4,heq2),fs,'./Pistas/pista_4_TOT_ls.wav');
	wavwrite(conv(aud_ea5,heq2),fs,'./Pistas/pista_5_TOT_ls.wav');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if(bool_ej(3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Ejercicio 2 %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diagrama de polos y ceros
figure
b = h;	% Coeficientes del numerador de H(z)
a = 1;	% Coeficientes del denominador de H(z)
zplane(b,a)



% Hallo los ceros de H.
ceros = roots(h);

% Destaco los que se encuentran afuera
aux = ceros.*(ceros>1);
aux2 = ceros-aux;

% Elimino los 0's redundantes
ceros_no_min =0;
ceros_min = 0;
for i=1:length(aux)
	if(aux(i))
		ceros_no_min = [ceros_no_min ; aux(i)];
	else
		ceros_min = [ceros_min ; aux2(i)];
	end
end
ceros_no_min = ceros_no_min(2:end);
ceros_min = ceros_min(2:end);

% Pongo los ceros dentro del círculo unitario
ceros_inv_conj = 1./conj(ceros_no_min);
fase_min = [ceros_min ; ceros_inv_conj];

h_fm = (poly(fase_min));

figure
zplane(1,h_fm);
xlabel('Real');
ylabel('Imaginario');


% Genero Hiir (compensador)
k = abs(h(1)*prod(ceros_no_min));
Hiir = 1./freqz(h_fm,1,nfft)/k;
Hiir = Hiir(1:end/2);

H = fft(h,nfft)';
H = H(1:end/2);

figure
zplane(h,h_fm);
xlabel('Real');
ylabel('Imaginario');

figure
plot(w,mag2db(abs((Hiir).*H)),'LineWidth',4,'r')
hold on
plot(w,mag2db(abs(H)),'LineWidth',4)
plot(w,mag2db(abs(Hiir)),'LineWidth',4,'color',[0 0.5 0])

lab_x = ["Frecuencia $\\frac{w}{\\pi}$"];
lab_y = ["Amplitud [\\si{\\dB}]"];
leyenda = ["Respuesta en frecuencia del sistema ecualizado"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema electroacústico"];
leyenda = [leyenda; "Respuesta en frecuencia del filtro ecualizador"];
loc = 'South';
AXIS = [0 1 -15 15];

xlabel(lab_x)
ylabel(lab_y)
grid minor

%set_graph('plot',[lab_x; lab_y],leyenda, loc, AXIS, 1);

figure
hold on
plot(w,grpdelay(1,h_fm,nfft/2),'LineWidth',4,'r');
plot(w,grpdelay(h,1,nfft/2),'LineWidth',4)
plot(w,grpdelay(h,h_fm,nfft/2),'LineWidth',4,'color',[0 .5 0]);

lab_y3 = ["Retardo de fase"];
leyenda = ["Filtro ecualizador"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema EA"];
leyenda = [leyenda; "Respuesta en frecuencia del sistema total"];
loc = 'North';
AXIS = [0 1 -20 080];

set_graph('plot',[lab_x; lab_y3],leyenda, loc, AXIS, 1);

heq3=filter(1,h_fm,[1 zeros(1,length(h)-1)]);

if(bool_audio)
	wavwrite(conv(aud_ea1,heq3),fs,'./Pistas/pista_1_TOT_iir.wav');
	wavwrite(conv(aud_ea2,heq3),fs,'./Pistas/pista_2_TOT_iir.wav');
	wavwrite(conv(aud_ea3,heq3),fs,'./Pistas/pista_3_TOT_iir.wav');
	wavwrite(conv(aud_ea4,heq3),fs,'./Pistas/pista_4_TOT_iir.wav');
	wavwrite(conv(aud_ea5,heq3),fs,'./Pistas/pista_5_TOT_iir.wav');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
