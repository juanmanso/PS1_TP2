clc;clear all; close all;
set(0,'DefaultFigureWindowStyle','docked')
%set(legend,'Interpreter','latex')
%format long;%format;%resetea el formato
%% TP Procesamiento Seniales 1
%           Estimacion LS
% FIUBA (Facultad de Ingenieria, Universidad de Buenos Aires)
%
% Realizado por
%   ·Eichenbaum, Daniel P.95233     Email: leinad_reverso@yahoo.com.ar
%   ·Accifonte, Franco  P.93799     Email: franco.accifonte@gmail.com
%---------------------------------------------------------- 
%% Variables globales
% Generacion del modelo
    load ('ensayo.mat')
    datos=datos(1:end,:);
    tita=tita(1:end);
    g=9.8;%m/s^2
    ax=datos(:,1);
    ay=datos(:,2);
    stita=sin(tita);
    ctita=cos(tita);
    N=length(tita);
% Generacion de la trayectoria
    load ('puntos.mat');
    pA=A;pB=B;pC=C;pD=D;
    load('acel.mat');

%% 1. Generacion del modelo
% Y=Ax+eta
% Ver el informe para entender de donde viene
Y=[ax+g*stita; ay+g*ctita];
A=[-g*stita ones(N,1) zeros(N,2);zeros(N,2) -g*ctita ones(N,1)];
x=zeros(4,1);%[s_x B_x s_y B_y]


%% 2. Estimacion de parametros
%x_hat=[     escalax;
%            sesgox;
%            escalay;
%            sesgox];

%% 2.1 Estimacion de la varianza del ruido
C_eta=[0.25 0;0 0.64];
G=inv(C_eta);

AHporG(:,[1 2])=A(:,[1 2]).*G(1,1);
AHporG(:,[3 4])=A(:,[3 4]).*G(2,2);
AHporG=AHporG';

x_hat=inv(AHporG*A)*AHporG*Y; %Con producto interno
x_hat_tradicional=inv(A'*A)*A'*Y;   %Sin nada
x_hat_moore=pinv(A)*Y;  % DVS (norma de X mínima)

%% 2.2 Varianza del estimador
C_xhat=inv( AHporG * A );

%% 3. Estimacion de la trayectoria
% El objetivo de este inciso 
% Será la de Estimar La posicion final.
% Para ello se debe integrar la aceleracion
% dos veces.

%% 3.1 Estimacion de la aceleracion real
    xAccReal=(Aerr(:,1)-x_hat(2))/(x_hat(1)+1);
    yAccReal=(Aerr(:,2)-x_hat(4))/(x_hat(3)+1);
    figure;
    subplot(2,1,1); %Acc x
        plot(t,xAccReal)
        ylabel('a_x[m/s^2]');xlabel('t');
        ylim([-2 2]);
    subplot(2,1,2); %Acc y
        plot(t,yAccReal);
        ylabel('a_y[m/s^2]');xlabel('t');
        ylim([-3 3]);
    figure;
        Ts=t(2)-t(1);
        Fs=1/Ts;
        winLength=0.5;
        wa=kaiser(fix(winLength/Ts),0.5);
    subplot(2,1,1);
        spectrogram(xAccReal,wa,fix(length(wa)*0.5),512,Fs,'yaxis')
        ylim([0 10]);
    subplot(2,1,2);
        spectrogram(yAccReal,wa,fix(length(wa)*0.5),512,Fs,'yaxis')
        ylim([0 10]);
        colormap(hsv)
%% 3.1b Filtrado de la aceleracion: Filtro el ruido con un FIR y aplico el mismo esquema que antes.
%    Aerr(:,1)=sin(2*pi*t*15);% Señal de testeo
    n = 100;
    f=10;
    ohm = [0 2*f*Ts 1]; %Se divide por pi porque esta normalizado
    a = [1 0];
    up = [1.01 0.01];
    lo = [0.99 -0.01];
    b = fircls(n,ohm,a,up,lo);%,'both');
    xAccFiltrada=filtfilt(b,1,Aerr(:,1)); %Filtro y desplazo M/2      %Aplico el filtro a las aceleraciones medidas
    yAccFiltrada=filtfilt(b,1,Aerr(:,2));
    
    figure;
    subplot(2,1,1); %Acc x
        plot(t,xAccFiltrada)
        ylabel('a_x[m/s^2]');xlabel('t');
        ylim([-2 2]);
    subplot(2,1,2); %Acc y
        plot(t,yAccFiltrada);
        ylabel('a_y[m/s^2]');xlabel('t');
        ylim([-3 3]);
    figure;
        Ts=t(2)-t(1);
        Fs=1/Ts;
        winLength=0.5;
        wa=kaiser(fix(winLength/Ts),0.5);
    subplot(2,1,1);
        spectrogram(xAccFiltrada,wa,fix(length(wa)*0.5),512,Fs,'yaxis')
        ylim([0 20]);
    subplot(2,1,2);
        spectrogram(yAccFiltrada,wa,fix(length(wa)*0.5),512,Fs,'yaxis')
        ylim([0 20]);
        colormap(hsv)

%% 3.2 Estimacion de la velocidad Real
xVelRaw=cumtrapz(t,Aerr(:,1))+1;
yVelRaw=cumtrapz(t,Aerr(:,2))+3;

xVelReal=cumtrapz(t,xAccReal)+1;
yVelReal=cumtrapz(t,yAccReal)+3;

xVelFiltrada=cumtrapz(t,xAccFiltrada)+1;
yVelFiltrada=cumtrapz(t,yAccFiltrada)+3;

    figure;
    subplot(2,1,1); %vel x
        plot(t,xVelRaw,'linewidth',2);hold on;
        plot(t,xVelReal,'linewidth',2);
        plot(t,xVelFiltrada,'linewidth',2);hold off;
        ylabel('v_x[m/s]');xlabel('t');legend({'Raw','$\hat{V}_{real}$','$\hat{V}_{filtrada}$'},'Interpreter','Latex');
        ylim([0 15]);
    subplot(2,1,2); %vel y
        plot(t,yVelRaw,'linewidth',2);hold on;
        plot(t,yVelReal,'linewidth',2);
        plot(t,yVelFiltrada,'linewidth',2);hold off;
        ylabel('v_y[m/s]');xlabel('t');
        ylim([-5 5]);
        
%% 3.3 Estimacion de la Posicion Real
xPosRaw=cumtrapz(t,xVelRaw);
yPosRaw=cumtrapz(t,yVelRaw);

xPosReal=cumtrapz(t,xVelReal);
yPosReal=cumtrapz(t,yVelReal);

xPosFiltrada=cumtrapz(t,xVelFiltrada);
yPosFiltrada=cumtrapz(t,yVelFiltrada);




    ha=figure;
        plot(xPosRaw,yPosRaw,'b','LineWidth',2);hold on;
        plot(xPosReal,yPosReal,'r','LineWidth',2);
        plot(xPosFiltrada,yPosFiltrada,'g','LineWidth',2);
                
    radios=[sqrt(C_xhat(2,2))+sqrt(C_xhat(1,1))*max(abs(xAccReal))%+sqrt(.25)   %Descomentar estas dos lineas para agregar la varianza del ruido blanco
            sqrt(C_xhat(4,4))+sqrt(C_xhat(3,3))*max(abs(yAccReal))];%+sqrt(.64)];
    radios=radios*t(end)*t(end)*3/2*1./(1+[x_hat(1);x_hat(3)]);
    
    redondo([xPosRaw(end) yPosRaw(end)],radios,ha,'--b');
    redondo([xPosReal(end) yPosReal(end)],radios,ha,'--r');
    redondo([xPosFiltrada(end) yPosFiltrada(end)],radios,ha,'--g');
    legend('Trayectoria del vehiculo sin correcciones'...
            ,'Trayectoria del vehiculo con correciones'...
            ,'Trayectoria con aceleracion filtrada'...
            ,'Intervalo de confianza para la posición final'...
            ,'Intervalo de confianza para la posición final'...
            ,'Intervalo de confianza para la posición final','Location','northwest');
    text(pA(1),pA(2),'A');
    text(pB(1),pB(2),'B');
    text(pC(1),pC(2),'C');
    text(pD(1),pD(2),'D');
    axis equal;
    grid on
    axis([-10 650 -50 240])

    
%% Punto 4
%% 4.1 Cantidad de puntos mimios

puntos=[pA;pB;pC;pD];
dists=0;
for ii=1:4
    for jj=ii:4
        if ii~=jj
            dists(end+1)=norm(puntos(ii,:)-puntos(jj,:));
        end
    end
end
dists=dists(2:end);


Nx=t(end)*t(end)*3/2;
Nx=Nx*sqrt(C_eta(1,1));
Nx=Nx/min(dists);
Nx=Nx*(1+max(abs(xAccReal)/9.8));
Nx=Nx*Nx;

Ny=t(end)*t(end)*3/2;
Ny=Ny*sqrt(C_eta(2,2));
Ny=Ny/min(dists);
Ny=Ny*(1+max(abs(yAccReal)/9.8));
Ny=Ny*Ny;