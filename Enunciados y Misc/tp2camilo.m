clear all;close all;

% Cargo las struct de los .mat %
Ensayo=load("Material/ensayo.mat");
Puntos=load("Material/puntos.mat");
Acel=load("Material/acel.mat");

% Paso los datos de struct a matrices %
  
  % Datos de aceleracion y angulo (Ensayo) %
  Tita=Ensayo.tita; % angulos en radianes %
  N=length(Tita);
  Amx=Ensayo.datos(:,1);
  Amy=Ensayo.datos(:,2);
  
  % Puntos posibles de llegada %
  PuntosLlegadaXY=zeros(4,2);
  PuntosLlegadaXY(1,:)=Puntos.A;
  PuntosLlegadaXY(2,:)=Puntos.B;
  PuntosLlegadaXY(3,:)=Puntos.C;
  PuntosLlegadaXY(4,:)=Puntos.D;

  % Aceleracion de trayectoria%
  Tiempo=Acel.t;
  T=length(Tiempo);
  AcelVehiculoX=Acel.Aerr(:,1);
  AcelVehiculoY=Acel.Aerr(:,2);
  
% Obtengo los datos iniciales y de ruido (datos.txt) %
N=length(Tita);
PosInicial=[0 0];
VelInicial=[1 3];
VarRuido=[0.25 0.64];
g=9.8; % modulo de la aceleracion de la gravedad en m/s2 %


%% Modelo en X: Amx=Gx*cx+Vx %

  % Armo la matriz Gx para el modelo en la direccion X %

  Arx=-g.*sin(Tita); % vector de valores de Arx %
  Gx=[ones(N,1) Arx];

  % Obtengo el cx que minimiza el ECM %

  cx=inv((Gx')*Gx)*(Gx')*Amx;

  % Obtengo los valores estimados de Esx y Eex %
  
  Esx=cx(1) % Error de sesgo en X %
  kx=cx(2);
  Eex=kx-1 % Error de escala en X % 

  % Grafico para visualizar el ajuste %
  figure;
  hold on;
  plot(Arx,Amx,'bo');

  Xfit=kx.*Arx+Esx;
  plot(Arx,Xfit,'r','Linewidth',3);
  title('Ajuste para mediciones en direccion X');
  ylabel('AcelX');
  xlabel('-|g|\cdotsen(\theta)');
  
  figure;
  plot(Tita,Amx);
  hold on;
  plot(Tita,Xfit,'r','Linewidth',3);
  title('Ajuste para mediciones en direccion X');
  ylabel('AcelX');
  xlabel('\theta (rad)');
  xlim([0 6.28])
  
  % Varianza del estimador cx %
  
  CovEstimadorX=VarRuido(1)*inv((Gx')*Gx)

%% Modelo en Y: Amy=Gy*cy+Vy
  
  % Armo la matriz Gy para el modelo en la direccion Y %

  Ary=-g.*cos(Tita); % vector de valores de Ary %
  Gy=[ones(N,1) Ary];

  % Obtengo el cy que minimiza el ECM %

  cy=inv((Gy')*Gy)*(Gy')*Amy;

  % Obtengo los valores estimados de Esy y Eey %
  
  Esy=cy(1) % Error de sesgo en Y %
  ky=cy(2);
  Eey=ky-1 % Error de escala en Y % 

  % Grafico para visualizar el ajuste %
  figure;
  hold on;
  plot(Ary,Amy,'bo');

  Yfit=ky.*Ary+Esy;
  plot(Ary,Yfit,'r','Linewidth',3);
  title('Ajuste para mediciones en direccion Y');
  ylabel('AcelY');
  xlabel('-|g|\cdotcos(\theta)');
  
  figure;
  plot(Tita,Amy);
  hold on;
  plot(Tita,Yfit,'r','Linewidth',3);
  title('Ajuste para mediciones en direccion Y');
  ylabel('AcelY');
  xlabel('\theta (rad)');
  xlim([0 6.28])
  
  % Varianza del estimador cy %
  
  CovEstimadorY=VarRuido(2)*inv((Gy')*Gy)
   
%% Calculo trayectorias (corregida y sin corregir) 
   
  % Corrijo la aceleracion en X e Y%
  
  AcelCorrX=(AcelVehiculoX.-Esx)./(1-Eex);
  AcelCorrY=(AcelVehiculoY.-Esy)./(1-Eey);

  % Calculo las trayectorias corregidas %

  VelX=VelInicial(1).+cumtrapz(Tiempo,AcelCorrX);
  PosX=PosInicial(1).+cumtrapz(Tiempo,VelX);

  VelY=VelInicial(2).+cumtrapz(Tiempo,AcelCorrY);
  PosY=PosInicial(2).+cumtrapz(Tiempo,VelY);
  
  figure;
  plot(PosX,PosY,'.');
  hold on;

  plot(PuntosLlegadaXY(1,1),PuntosLlegadaXY(1,2),'rx','Linewidth',2);
  plot(PuntosLlegadaXY(2,1),PuntosLlegadaXY(2,2),'mx','Linewidth',2);
  plot(PuntosLlegadaXY(3,1),PuntosLlegadaXY(3,2),'kx','Linewidth',2);
  plot(PuntosLlegadaXY(4,1),PuntosLlegadaXY(4,2),'gx','Linewidth',2);
  title('Trayectoria corregida');
  xlabel('Posicion eje X');
  ylabel('Posicion eje Y');
  legend('Trayectoria Corregida','Punto A','Punto B','Punto C','Punto D',"location","northwest");
  
  Tf=Tiempo(T);
  dX=Tf*Tf*0.5*3*(sqrt(CovEstimadorX(1,1))+max(abs(AcelCorrX))*sqrt(CovEstimadorX(2,2)));
  dY=Tf*Tf*0.5*3*(sqrt(CovEstimadorY(1,1))+max(abs(AcelCorrY))*sqrt(CovEstimadorY(2,2)));
  for i=1:4
    elipse(dX,dY,PuntosLlegadaXY(i,1),PuntosLlegadaXY(i,2),i)
  end
  % Calculo las trayectorias sin corregir para comparar %

  VelXErr=VelInicial(1).+cumtrapz(Tiempo,AcelVehiculoX);
  PosXErr=PosInicial(1).+cumtrapz(Tiempo,VelXErr);

  VelYErr=VelInicial(2).+cumtrapz(Tiempo,AcelVehiculoY);
  PosYErr=PosInicial(2).+cumtrapz(Tiempo,VelYErr);

  figure;
  plot(PosXErr,PosYErr,'.');
  hold on;
  plot(PuntosLlegadaXY(1,1),PuntosLlegadaXY(1,2),'rx','Linewidth',2);
  plot(PuntosLlegadaXY(2,1),PuntosLlegadaXY(2,2),'mx','Linewidth',2);
  plot(PuntosLlegadaXY(3,1),PuntosLlegadaXY(3,2),'kx','Linewidth',2);
  plot(PuntosLlegadaXY(4,1),PuntosLlegadaXY(4,2),'gx','Linewidth',2);

  title('Trayectoria sin corregir');
  xlabel('Posicion eje X');
  ylabel('Posicion eje Y');
  legend('Trayectoria sin corregir','Punto A','Punto B','Punto C','Punto D',"location","northwest");
  
  
  
  
