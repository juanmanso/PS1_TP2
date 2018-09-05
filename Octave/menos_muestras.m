
config_m;

%%%%%%%% Trabajo Práctico 2 %%%%%%%%%

N=2222;

% Se cargan las simulaciones
ensayo_struct = load('ensayo.mat');
puntos_struct = load('puntos.mat');
acel_struct = load('acel.mat');

% Datos adicionales de la simulación
p_inicial = [0 0];
v_inicial = [1 3];
sigma_ruido = [0.25 0.64];


%% Variables propias
g = 9.8;



%%%%%% ESTIMACIÓN DEL SESGO Y FACTOR DE ESCALA %%%%%%%
theta = ensayo_struct.tita;
a_med = ensayo_struct.datos;
l_viejo = length(theta);
l = floor(l_viejo/N);
theta=theta(1:l:end);
a_med=a_med(1:l:end,:);
l = length(theta);


Y_x = [a_med(:,1)+g.*sin(theta)];
Y_y = [a_med(:,2)+g.*cos(theta)];
A_x = [ones(l,1) -g.*sin(theta)];
A_y = [ones(l,1) -g.*cos(theta)];

A = [A_x A_x*0; A_y*0 A_y];
Y = [Y_x; Y_y];

X_x = ((A_x'*A_x)^(-1))*(A_x')*Y_x;
X_y = ((A_y'*A_y)^(-1))*(A_y')*Y_y;

bx = X_x(1);
sx = X_x(2);
by = X_y(1);
sy = X_y(2);

covx = sigma_ruido(1) * (A_x'*A_x)^(-1);
covy = sigma_ruido(2) * (A_y'*A_y)^(-1);
sigma_bx = covx(1,1);
sigma_sx = covx(2,2);
sigma_by = covy(1,1);
sigma_sy = covy(2,2);

sigma_b = [sigma_bx, sigma_by];
sigma_s = [sigma_sx, sigma_sy];

sigmas = [sigma_b, sigma_s, sigma_ruido];

%%%%%%% CÁLCULO DE LA TRAYECTORIA %%%%%

t = acel_struct.t;
a_tray = acel_struct.Aerr;

PA = puntos_struct.A;
PB = puntos_struct.B;
PC = puntos_struct.C;
PD = puntos_struct.D;

P = [PA;PB;PC;PD];

errores = [a_tray(:,1)*sx+bx a_tray(:,2)*sy+by];

v_calc = v_inicial + cumtrapz(t,a_tray) - cumtrapz(t,errores);
p_calc = p_inicial + cumtrapz(t,v_calc) - cumtrapz(t,errores);

p_sin_cor = cumtrapz(t,(cumtrapz(t,a_tray)+v_inicial))+p_inicial;

figure
plot(p_calc(:,1),p_calc(:,2),'LineWidth',2)
hold on
plot(p_sin_cor(:,1),p_sin_cor(:,2),'LineWidth',2,'color',myGreen);

deltap=3*(t.^2/2.*(sigma_b+max(a_tray).*sigma_s)+sigma_ruido);
plot(p_calc(:,1)+deltap(:,1),p_calc(:,2)+deltap(:,2),'c.');
plot(p_calc(:,1)-deltap(:,1),p_calc(:,2)-deltap(:,2),'c.');
plot(P(:,1),P(:,2),'rx','Markersize',03)

xlab='Coordenada X de posición';
ylab='Coordenada Y de posición';
%graficar_incertidumbre(P,t,a_tray,sigmas)
leyenda=['Trayectoria realizada por el vehículo';'Trayectoria sin correción';'Incerteza de la estimación'];
loc='Northeast';
ejes=[0 600];

set_graph('plot',[xlab;ylab],leyenda,loc,ejes,1);



