#! /usr/bin/octave -qf


#**************************************
# Script de Octave para calcular las tensiones y corrientes del boost en base al Duty aplicado y otros parÁmetros
#******************
# Los parametros que recibe son booleanos para -p (imprimir figuras) y -results (exportar resultados)
# 
# Uso:
#	$ ./graph -p <bool> -results <bool>
#
#	_____________________________________________________________
#	COMPROBAR LA VALIDACIÓN DE ARGUMENTOS EN LA CARPETA FUNCIONES
#	
#**************************************

config_graph

#******************
# Validacion de argumentos
#******************

if(!validate_arguments(argv(),nargin))
	exit(1);
endif
arg_list = argv();

#******************
# Macros
#******************
DELIMITER='|';
ARG_PRINT_POS=2;
ARG_RESULTS_POS=4;
POSICION_MED=1;		% Posición en el arreglo de formato del grafico 
POSICION_AJUSTE=2;	% correspondiente a las mediciones y las curvas de ajuste
MAX_TRANSISTOR=3;
AMPERE_UNITS=1;
mAMPERE_UNITS=1E3;
CURRENT_FACTOR=mAMPERE_UNITS;	% Para calcular y graficar en A ==> 1; mA ==> 1e3


#******************
# Variables de opciones | Variables booleanas
#******************
bool_print=str2num(arg_list{ARG_PRINT_POS});		% Si vale 1, se imprimen los gráficos
bool_export_results=str2num(arg_list{ARG_RESULTS_POS});	% Si vale 1, genera un archivo con los resultados (arreglos)
bool_close_all=0;	% Si vale 1, al finalizar el script, cierra las figuras
bool_grid_minor=1;
bool_sat_mA_units=1;	% Si vale 1, las curvas de salida se grafican en mA (independientemente del CURRENT_FACTOR
bool_bc337_40=0;	% Si vale 1, se usa la simulación con el tr BC337-40. Si vale 0, se usa el BC337

print_color="-color";	# Si quiere imprimirse color, -color. Si no, -mono
format_out_name=".tex";			% Configuro el formato de impresión
format_print_config="-depslatex";	% Idem
print_path="./";		% Ruta de salida


#******************
# Nombres de los archivos
#******************
med_Ic_Vbe_tr1_file="./Datos/Mediciones/t1_______ic(vbe).txt";
med_Ic_Vbe_tr2_file="./Datos/Mediciones/t2_______ic(vbe).txt";
med_Ic_Vbe_tr3_file="./Datos/Mediciones/t3_______ic(vbe).txt";
med_Ic_VcE_7_5mA_tr1_file="./Datos/Mediciones/t1_ic(vce)_7_5mA.txt";
med_Ic_VcE_7_5mA_tr2_file="./Datos/Mediciones/t2_ic(vce)_7_5mA.txt";
med_Ic_VcE_7_5mA_tr3_file="./Datos/Mediciones/t3_ic(vce)_7_5mA.txt";
med_Ic_VcE_15mA_tr1_file="./Datos/Mediciones/t1_ic(vce)__15mA.txt";
med_Ic_VcE_15mA_tr2_file="./Datos/Mediciones/t2_ic(vce)__15mA.txt";
med_Ic_VcE_15mA_tr3_file="./Datos/Mediciones/t3_ic(vce)__15mA.txt";

sim_IC_VBE_file="./Datos/Simulaciones/BC337/SIM_______Ic(Vbe).txt";
sim_IC_VCE_7_5mA_file="./Datos/Simulaciones/BC337/SIM_Ic(Vce)_7_5mA.txt";
sim_IC_VCE_15mA_file="./Datos/Simulaciones/BC337/SIM_Ic(Vce)__15mA.txt";

sim_nuestro_IC_VBE_file="./Datos/Modelo_propio/NUESTRO_MODEL_______Ic(Vbe).txt";
sim_nuestro_IC_VCE_7_5mA_file="./Datos/Modelo_propio/NUESTRO_MODEL_Ic(Vce)_7_5mA.txt";
sim_nuestro_IC_VCE_15mA_file="./Datos/Modelo_propio/NUESTRO_MODEL_Ic(Vce)__15mA.txt";


#******************
# Carga de las simulaciones (fabricante y modelo propio)
#******************
SIM_IC_VS_VBE=dlmread(sim_IC_VBE_file,'\t',1,0); %(V [V], I[A])
	SIM_IC_VS_VBE(:,2)=SIM_IC_VS_VBE(:,2).*CURRENT_FACTOR;
SIM_ICSAT_7_5mA=dlmread(sim_IC_VCE_7_5mA_file,'\t',1,0);  %(V [V], I[A])
%	SIM_ICSAT_7_5mA(:,2)=SIM_ICSAT_7_5mA(:,2).*CURRENT_FACTOR;
SIM_ICSAT_15mA=dlmread(sim_IC_VCE_15mA_file,'\t',1,0);  %(V [V], I[A])
%	SIM_ICSAT_15mA(:,2)=SIM_ICSAT_15mA(:,2).*CURRENT_FACTOR;

SIM_nuestro_IC_VS_VBE=dlmread(sim_nuestro_IC_VBE_file,'\t',1,0); %(V [V], I[A])
	SIM_nuestro_IC_VS_VBE(:,2)=SIM_nuestro_IC_VS_VBE(:,2).*CURRENT_FACTOR;
SIM_nuestro_ICSAT_7_5mA=dlmread(sim_nuestro_IC_VCE_7_5mA_file,'\t',1,0); %(V [V], I[A])
%	SIM_nuestro_ICSAT_7_5mA(:,2)=SIM_nuestro_ICSAT_7_5mA(:,2).*CURRENT_FACTOR;
SIM_nuestro_ICSAT_15mA=dlmread(sim_nuestro_IC_VCE_15mA_file,'\t',1,0); %(V [V], I[A])
%	SIM_nuestro_ICSAT_15mA(:,2)=SIM_nuestro_ICSAT_15mA(:,2).*CURRENT_FACTOR;


files=[	med_Ic_Vbe_tr1_file; med_Ic_Vbe_tr2_file; med_Ic_Vbe_tr3_file;
	med_Ic_VcE_7_5mA_tr1_file; med_Ic_VcE_7_5mA_tr2_file; med_Ic_VcE_7_5mA_tr3_file;
	med_Ic_VcE_15mA_tr1_file; med_Ic_VcE_15mA_tr2_file; med_Ic_VcE_15mA_tr3_file];


#******************
# Arreglos que almacenarán los parámetros calculados
#******************
Is_array=zeros(1,MAX_TRANSISTOR);
Vth_array=zeros(1,MAX_TRANSISTOR);
ro_array=zeros(2*MAX_TRANSISTOR,2); %La 1er columna es ro y la segunda lambda, cada fila es un tr, cada 3 filas cambia el ICsat 
VA_array=zeros(2,MAX_TRANSISTOR);   % La primer fila es para ICsat=5mA, la segunda, ICsat=20mA
BETA_array=[280,286,289];	%Beta's medidos con el multímetro p/cada tr
BETA_fab=292;
BETA_nuestra_sim=289;  


#******************
#Cadenas de caracteres para los gráficos
#******************
ley_med_tr1="Curva experimental $tr_1$";
ley_med_tr2="Curva experimental $tr_2$";
ley_med_tr3="Curva experimental $tr_3$";
ley_aj_tr1="";
ley_aj_tr2="";
ley_aj_tr3="";
ley_sim="Curva fabricante";
ley_sim_nuestra="Curva parametrizada";
transf_title="";%"Curva de transferencia";
tr_out_20mA_title="";%"Curva de salida 20mA";
tr_out_5mA_title="";%"Curva de salida 5mA";
transconduct_title="";%"Curva de transconductancia";
rpi_title="";%"Curva de R_{PI}";
output_file_name1="graf_ic(vbe)";
output_file_name2="graf_Icsat=7_5mA";
output_file_name3="graf_Icsat=15mA";
output_file_name4="graph_gm";
output_file_name5="graph_rpi";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						MAIN							%				
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ajuste_subindice=[8,7,12,	13,15,16,	16,16,24]; #en vez de 24 podria ir 19
superior_subindice=[0,0,0,	0,0,0,	0,0,0];
tr=1;				% Inicializo los parámetros de las iteraciones
med=1;
file_aux=(tr+3*(med-1));	% file_aux es la posición de los archivos para cada tr de cada medicion

figure(1)
hold on

if(bool_grid_minor)
grid minor
else
grid on
endif

	while(tr<=MAX_TRANSISTOR)  % Se calcula Is y Vth para los 3 transistores

		k=ajuste_subindice(file_aux);
		f=superior_subindice(file_aux);
		BJT=dlmread(files(file_aux,:),DELIMITER,1,0); %(V [V], I[mA])
		BJT(2,:)=BJT(2,:).*1e-3;		  %(V [V], I[A])

		lBJT=BJT;				
		lBJT(2,:)=log(lBJT(2,:)); 		% ln(ic)=ln(IS)+vbe/Vth
% A = 1/Vth
% B = ln(Is)
% Condiciones iniciales del ajuste elegidas: Vth=35mV y Is=1e-13A
% => A = 1/25e-3 = 40
% => B = ln(1e-13) = -29.9
		Param=fmins('A_x_mas_B',[1/30e-3 log(1e-13)],[0,1e-3,0,0,0,0],[],lBJT(:,k:end-f));
		Is_array(tr)=exp(Param(2));           	%Is=exp(ln(Is)) [A]
		Vth_array(tr)=1/Param(1);		%Vth=1/(1/Vth)  [V] ~35mV

%BJT(1,k)
%BJT(1,end)
			############################################################
		% Para graficar las curvas TEORICAS, nos ayudamos de un vector auxiliar para las tensiones
		V=BJT(1,k):0.001:BJT(1,end-f);		


	#Configuro el gráfico y calculo gm cuya longitud depende de las muestras
	#La primer fila del gm_tr? tiene los valores de gm[A/V] y la segunda Ic[A]
		if(tr==1)				
			gm_tr1=zeros(2,length(BJT)-1);
			gm_tr1(1,:)=diff(BJT(2,:))./diff(BJT(1,:));
			gm_tr1(2,:)=BJT(2,2:end);
			graph_line_format=["bo";"b-"];
		elseif(tr==2)
			gm_tr2=zeros(2,length(BJT)-1);
			gm_tr2(1,:)=diff(BJT(2,:))./diff(BJT(1,:));
			gm_tr2(2,:)=BJT(2,2:end);
			graph_line_format=["ro";"r-"];
		else
			gm_tr3=zeros(2,length(BJT)-1);
			gm_tr3(1,:)=diff(BJT(2,:))./diff(BJT(1,:));
			gm_tr3(2,:)=BJT(2,2:end);
			graph_line_format=["go";"g-"];
		endif

	#Grafico primero las mediciones y luego los ajustes
		semilogy(BJT(1,1:end),BJT(2,1:end).*CURRENT_FACTOR,graph_line_format(POSICION_MED,:),'Markersize',5)
		semilogy(V(:),CURRENT_FACTOR.*exp(V(:)./Vth_array(tr)).*(Is_array(tr)),graph_line_format(POSICION_AJUSTE,:),'Linewidth',2)
		
			############################################################

		tr++;
		file_aux=(tr+3*(med-1)); %Es equivalente a file_aux++, pero se pierde el sentido de lo que significa
	endwhile

%Grafico de simulacion de PHIL_BJT.lib
	semilogy(SIM_IC_VS_VBE(:,1),(SIM_IC_VS_VBE(:,2)),'m-','Linewidth',7)
%Grafico de simulacion de nuestro tr
	semilogy(SIM_nuestro_IC_VS_VBE(1:5:end,1),(SIM_nuestro_IC_VS_VBE(1:5:end,2)),'y-','Linewidth',7)	
hold off
axis([0.56 0.750 10e-5*CURRENT_FACTOR .1*CURRENT_FACTOR]);

xlabel("$V_{BE}$ [V]")

	if(CURRENT_FACTOR==AMPERE_UNITS)
		ylabel("$I_{C}$ [A]")
	else
		ylabel("$I_{C}$ [mA]")
	endif

s_param1=sprintf(" | $I_S=\\SI{%2.3g}{\\A}$ - $V_{th}=\\SI{%4.3g}{\\mV}$",Is_array(1),Vth_array(1)*1e3);
s_param2=sprintf(" | $I_S=\\SI{%2.3g}{\\A}$ - $V_{th}=\\SI{%4.3g}{\\mV}$",Is_array(2),Vth_array(2)*1e3);
s_param3=sprintf(" | $I_S=\\SI{%2.3g}{\\A}$ - $V_{th}=\\SI{%4.3g}{\\mV}$",Is_array(3),Vth_array(3)*1e3);
s_param4=sprintf(" | $I_S=\\SI{%2.3g}{\\A}$ - $V_{th}=\\SI{%4.3g}{\\mV}$",1.41e-13,27.4);
s_param5=sprintf(" | $I_S=\\SI{%2.3g}{\\A}$ - $V_{th}=\\SI{%4.3g}{\\mV}$",1.592e-13,26.78);
Hleg = legend(gca,[ley_med_tr1 s_param1],ley_aj_tr1,[ley_med_tr2 s_param2],ley_aj_tr2,[ley_med_tr3 s_param3],ley_aj_tr3,[ley_sim s_param4],[ley_sim_nuestra s_param5],"location","southeast");
legend('boxon');
set(Hleg,'FontName','Arial','FontSize',7);


title(transf_title)

if(bool_print)
print([print_path,output_file_name1,format_out_name],format_print_config, "-color");
endif



if(1)

%%%%%%%%%%%%%% Grafico de gm para cada tr en funcion de Id %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------------------------------------
% Grafico de la transconductancia gm calculada;		ln(ic)=ln(Is)+vbe/Vth
%  - A partir del cociente de diferencias
%  - A partir de gm=dic/dvbe=Ic/Vth
%---------------------------------------------------------------------------------------------------------------
I=0:0.1:(0.50*CURRENT_FACTOR);
gm_tr1=gm_tr1.*CURRENT_FACTOR;
gm_tr2=gm_tr2.*CURRENT_FACTOR;
gm_tr3=gm_tr3.*CURRENT_FACTOR;

%Cálculo gm de simulacion de PHIL_BJT.lib
gm_sim=diff(SIM_IC_VS_VBE(:,2))./diff((SIM_IC_VS_VBE(:,1)));
IC_gm_sim=SIM_IC_VS_VBE(2:end,2);

%Cálculo gm de simulacion de nuestro tr
gm_sim_nuestro=diff(SIM_nuestro_IC_VS_VBE(:,2))./diff((SIM_nuestro_IC_VS_VBE(:,1)));
IC_gm_sim_nuestro=SIM_nuestro_IC_VS_VBE(2:end,2);

figure(4)
  hold on
if(bool_grid_minor)
grid minor
else
grid on
endif
  plot(gm_tr1(2,:),gm_tr1(1,:),'bo','Markersize',5)
  plot(I,(I./Vth_array(1)),'b-','Linewidth',2)
  plot(gm_tr2(2,:),gm_tr2(1,:),'ro','Markersize',5)
  plot(I,(I./Vth_array(2)),'r-','Linewidth',2)
  plot(gm_tr3(2,:),gm_tr3(1,:),'go','Markersize',5)
  plot(I,(I./Vth_array(3)),'g-','Linewidth',2)

%Gráfico de simulacion de CD4007.lib
 
  plot(IC_gm_sim(1:5:end),gm_sim(1:5:end),'m-','Linewidth',7)

%Gráfico de simulación de nuestro tr

	  plot(IC_gm_sim_nuestro(1:5:end),gm_sim_nuestro(1:5:end),'y-','Linewidth',7)
	  
  hold off
if(CURRENT_FACTOR==AMPERE_UNITS)
	xlabel("$I_{C}$ [A]");
	ylabel("$g_m$ [S]");
else
	xlabel("$I_{C}$ [mA]");
	ylabel("$g_m$ [mS]");
endif

  s_param1=sprintf(" | $V_{th}=\\SI{%2.3g}{\\mV}$",Vth_array(1)*1e3);
  s_param2=sprintf(" | $V_{th}=\\SI{%2.3g}{\\mV}$",Vth_array(2)*1e3);
  s_param3=sprintf(" | $V_{th}=\\SI{%2.3g}{\\mV}$",Vth_array(3)*1e3);
  s_param4=sprintf(" | $V_{th}=\\SI{%4.3g}{\\mV}$",27.4);
  s_param5=sprintf(" | $V_{th}=\\SI{%4.3g}{\\mV}$",26.78);
  Hleg = legend(gca,[ley_med_tr1 s_param1],ley_aj_tr1,[ley_med_tr2 s_param2],ley_aj_tr2,[ley_med_tr3 s_param3],ley_aj_tr3,[ley_sim s_param4],[ley_sim_nuestra s_param5],"location","northwest");
  legend('boxon');
  set(Hleg,'FontName','Arial','FontSize',8);
  title(transconduct_title);
  axis([0,0.05,0,2.500].*CURRENT_FACTOR);

  if(bool_print)
  print([print_path,output_file_name4,format_out_name],format_print_config, "-color");
  endif
endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Grafico de Rpi para cada tr en funcion de Id %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------------------------------------
% Grafico de la transconductancia Rpi calculada;
%  - A partir del cálculo de gm
%---------------------------------------------------------------------------------------------------------------
I=0:0.1:(0.50*CURRENT_FACTOR);


figure(5)
  hold on
if(bool_grid_minor)
grid minor
else
grid on
endif
  plot(gm_tr1(2,:),BETA_array(1)./gm_tr1(1,:),'bo','Markersize',5)
  plot(I,BETA_array(1)./(I./Vth_array(1)),'b-','Linewidth',2)
  plot(gm_tr2(2,:),BETA_array(2)./gm_tr2(1,:),'ro','Markersize',5)
  plot(I,BETA_array(2)./(I./Vth_array(2)),'r-','Linewidth',2)
  plot(gm_tr3(2,:),BETA_array(3)./gm_tr3(1,:),'go','Markersize',5)
  plot(I,BETA_array(3)./(I./Vth_array(3)),'g-','Linewidth',2)

%Gráfico de simulacion de PHIL_BJT.lib
 
	plot(IC_gm_sim(30:end),BETA_fab./gm_sim(30:end),'m-','Linewidth',7)

%Gráfico de simulación de nuestro tr

	plot(IC_gm_sim_nuestro(1:end),BETA_nuestra_sim./gm_sim_nuestro(1:end),'y-','Linewidth',7)
	  
  hold off
if(CURRENT_FACTOR==AMPERE_UNITS)
	xlabel("$I_{C}$ [A]");
	ylabel("$R_{\\pi}$ [$\\SI{}{\\ohm}$]");
else
	xlabel("$I_{C}$ [mA]");
	ylabel("$R_{\\pi}$ [$\\SI{}{\\kilo\\ohm}$]");
endif

  s_param1=sprintf(" | $V_{th}=\\SI{%2.3g}{\\mV}$ - $\\beta$=%2.3g",Vth_array(1)*1e3,BETA_array(1));
  s_param2=sprintf(" | $V_{th}=\\SI{%2.3g}{\\mV}$ - $\\beta$=%2.3g",Vth_array(2)*1e3,BETA_array(2));
  s_param3=sprintf(" | $V_{th}=\\SI{%2.3g}{\\mV}$ - $\\beta$=%2.3g",Vth_array(3)*1e3,BETA_array(3));
  s_param4=sprintf(" | $V_{th}=\\SI{%4.3g}{\\mV}$ - $\\beta$=%2.3g",27.4,BETA_fab);
  s_param5=sprintf(" | $V_{th}=\\SI{%4.3g}{\\mV}$ - $\\beta$=%2.3g",26.76,BETA_nuestra_sim);
  Hleg = legend(gca,[ley_med_tr1 s_param1],ley_aj_tr1,[ley_med_tr2 s_param2],ley_aj_tr2,[ley_med_tr3 s_param3],ley_aj_tr3,[ley_sim s_param4],[ley_sim_nuestra s_param5],"location","northwest");
  legend('boxon');
  set(Hleg,'FontName','Arial','FontSize',8);
  title(rpi_title);
  axis([1E-3,50E-3, 0,10E-3].*CURRENT_FACTOR);

  if(bool_print)
  print([print_path,output_file_name5,format_out_name],format_print_config, "-color");
  endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

med++;			%Inicializo nuevamente los parámetros
tr=1;
file_aux=(tr+3*(med-1));
ro_aux=file_aux-3;

figure(2)
hold on

if(bool_grid_minor)
grid minor
else
grid on
endif

%Saco los valores de ro para tr1, tr2 y tr3 para 1mA y grafico	
	while(tr<=MAX_TRANSISTOR)

		k=ajuste_subindice(file_aux);
		BJT=dlmread(files(file_aux,:),DELIMITER,1,0); %(V [V], I[mA])	
		%BJT_1=fliplr(BJT_1);		  
		BJT_t=BJT';
		if(bool_sat_mA_units)
			CURRENT_FACTOR=mAMPERE_UNITS;
		endif
		
%-----------------------------------------------------      
% Ajuste de la corriente en la región de MAD
% El objetivo del ajuste mediante 'fmins' y la funcion auxiliar 'A_x_mas_B'
% es minimizar iterativamente el error cuadrático entre los valores de Ic 
% medidos y el resultado de A*VCE+B , donde A y B son los parametros de ajuste

% Condiciones iniciales del ajuste elegidas: A=1/ro=1/10 [kOhm] y B=Icsat=4.7 [mA]
% NOTA: observar que en el ajuste se utilizan ambas columnas, pero sólo se utilizan 
% algunas de las filas de los datos
%-----------------------------------------------------  
		Param=fmins('A_x_mas_B',[1/10E-3 4.7E-3],[0,1e-3,0,0,0,0],[],BJT_t(k:end,:)');	
				
		% Obtengo ro y lambda
		ro_array(ro_aux,1)=1/(Param(1)*1E-3);		%Ro=1/1/Ro
		ro_array(ro_aux,2)=Param(2)*1E-3;		%B=Icsat
		VA_array(1,tr)=Param(2)/Param(1);		% abs(VA)=Icsat*Ro
							
		V=BJT_t(k,1):0.0001:BJT_t(end,1);
		BJT_t(:,2)=BJT_t(:,2).*1e-3;

		if(tr==1)
			graph_line_format=["bo";"b-"];
		elseif(tr==2)
			graph_line_format=["ro";"r-"];
		else
			graph_line_format=["go";"g-"];
		endif

		plot(BJT_t(1:end,1),BJT_t(1:end,2).*CURRENT_FACTOR,graph_line_format(POSICION_MED,:),'Markersize',5)
		plot(V(:),((V(:)./ro_array(ro_aux,1))+ro_array(ro_aux,2)).*CURRENT_FACTOR,graph_line_format(POSICION_AJUSTE,:),'Linewidth',2)
							
		tr++;					
		file_aux=(tr+3*(med-1));
		ro_aux=file_aux-3;
	endwhile
%Grafico de simulacion de CD4007.lib
	plot(SIM_ICSAT_7_5mA(:,1),SIM_ICSAT_7_5mA(:,2).*CURRENT_FACTOR,'m-','Linewidth',7)
%Grafico de simulacion de nuestro tr
	plot(SIM_nuestro_ICSAT_7_5mA(:,1),SIM_nuestro_ICSAT_7_5mA(:,2).*CURRENT_FACTOR,'y-','Linewidth',7)
	
hold off	

xlabel("$V_{CE}$ [V]")

	if(CURRENT_FACTOR==AMPERE_UNITS)
		ylabel("$I_{C}$ [A]")
	else
		ylabel("$I_{C}$ [mA]")
	endif	

s_param1=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",ro_array(1+med-2,1)/1E3,ro_array(1+med-2,2)*1E3,VA_array(1,1));
s_param2=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",ro_array(2+med-2,1)/1E3,ro_array(2+med-2,2)*1E3,VA_array(1,2));
s_param3=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%4.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",ro_array(3+med-2,1)/1E3,ro_array(3+med-2,2)*1E3,VA_array(1,3));
s_param4=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",5.64,7.23,40.8);
s_param5=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",5.36,7.937,42.559);%ro_array(3+med-2,1)/1E3,ro_array(3+med-2,2)*1E3,VA_array(1,3));
Hleg = legend(gca,[ley_med_tr1 s_param1],ley_aj_tr1,[ley_med_tr2 s_param2],ley_aj_tr2,[ley_med_tr3 s_param3],ley_aj_tr3,[ley_sim s_param4],[ley_sim_nuestra s_param5],"location","southeast");
legend('boxon');
set(Hleg,'FontName','Arial','FontSize',8);
title(tr_out_5mA_title)
axis([0 3 0 9]);

if(bool_print)
print([print_path,output_file_name2,format_out_name],format_print_config, "-color");
endif
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

med++;			%Inicializo nuevamente los parámetros
tr=1;
file_aux=tr+3*(med-1);
ro_aux=file_aux-3;

figure (3)
hold on

if(bool_grid_minor)
grid minor
else
grid on
endif

	while(tr<=MAX_TRANSISTOR)

		k=ajuste_subindice(file_aux);
		BJT=dlmread(files(file_aux,:),DELIMITER,1,0);	
		%BJT_1=fliplr(BJT_1);							
		BJT(2,:)=BJT(2,:).*1e-3;		  %(V [V], I[A])
		BJT_t=BJT';
		if(bool_sat_mA_units)
			CURRENT_FACTOR=mAMPERE_UNITS;
		endif
%-----------------------------------------------------      
% Ajuste de la corriente en la región de MAD
% El objetivo del ajuste mediante 'fmins' y la funcion auxiliar 'A_x_mas_B'
% es minimizar iterativamente el error cuadrático entre los valores de Ic 
% medidos y el resultado de A*VCE+B , donde A y B son los parametros de ajuste

% Condiciones iniciales del ajuste elegidas: A=1/ro=1/1.5 [kOhm] y B=Icsat=17 [mA]
% NOTA: observar que en el ajuste se utilizan ambas columnas, pero sólo se utilizan 
% algunas de las filas de los datos
%-----------------------------------------------------  
		Param=fmins('A_x_mas_B',[1/1.5E-3 17E-3],[0,1e-3,0,0,0,0],[],BJT_t(k:end,:)');	
							
		% Obtengo ro y lambda
		ro_array(ro_aux,1)=1/Param(1);			%Ro
		ro_array(ro_aux,2)=(Param(2));			%Icsat
		VA_array(2,tr)=Param(2)/Param(1);		% abs(VA)=Icsat*Ro


		V=BJT_t(k,1):0.0001:BJT_t(end,1);

		if(tr==1)
			graph_line_format=["bo";"b-"];
		elseif(tr==2)
			graph_line_format=["ro";"r-"];
		else
			graph_line_format=["go";"g-"];
		endif

		plot(BJT_t(1:end,1),BJT_t(1:end,2).*CURRENT_FACTOR,graph_line_format(POSICION_MED,:),'Markersize',5)
		plot(V(:),((V(:)./ro_array(ro_aux,1))+ro_array(ro_aux,2)).*CURRENT_FACTOR,graph_line_format(POSICION_AJUSTE,:),'Linewidth',2)

												
		tr++;
		file_aux=tr+3*(med-1);
		ro_aux=file_aux-3;
	endwhile		
%Grafico de simulacion de CD4007.lib
	plot(SIM_ICSAT_15mA(:,1),SIM_ICSAT_15mA(:,2).*CURRENT_FACTOR,'m-','Linewidth',7)
%Grafico de simulacion de nuestro tr
	plot(SIM_nuestro_ICSAT_15mA(:,1),SIM_nuestro_ICSAT_15mA(:,2).*CURRENT_FACTOR,'y-','Linewidth',7)
hold off		
xlabel("$V_{CE}$ [V]")
	if(CURRENT_FACTOR==AMPERE_UNITS)
		ylabel("$I_{C}$ [A]")
	else
		ylabel("$I_{C}$ [mA]")
	endif
s_param1=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",ro_array(4,1)/1E3,ro_array(4,2)*1E3,VA_array(2,1));
s_param2=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",ro_array(5,1)/1E3,ro_array(5,2)*1E3,VA_array(2,2));
s_param3=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",ro_array(6,1)/1E3,ro_array(6,2)*1E3,VA_array(2,3));
s_param4=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",4.08,14.604,59.63);
s_param5=sprintf(" | $r_o=\\SI{%0.3g}{\\kilo\\ohm}$ - $I_{Cmad}$=\\SI{%0.3g}{\\mA} - $V_A=\\SI{%0.3g}{\V}$",1.65,15.073,24.9);%ro_array(6,1)/1E3,ro_array(6,2)*1E3,VA_array(2,3));
Hleg = legend(gca,[ley_med_tr1 s_param1],ley_aj_tr1,[ley_med_tr2 s_param2],ley_aj_tr2,[ley_med_tr3 s_param3],ley_aj_tr3,[ley_sim s_param4],[ley_sim_nuestra s_param5],"location","southeast");
legend('boxon');
set(Hleg,'FontName','Arial','FontSize',8);
title(tr_out_20mA_title)
axis([0 3 0 20]);

if(bool_print)
print([print_path,output_file_name3,format_out_name],format_print_config, "-color");
endif

if(bool_close_all)
close all;
endif

% Exportar resultados
if(bool_export_results)
dlmwrite("resultados_octave.txt","Resultados del análisis con Octave\nFila 1°: Arreglo Vth\t2°: Arreglo Is\t3°: Arreglo r0\t4°: Arreglo de ICsat","")
dlmwrite("resultados_octave.txt",[Vth_array;Is_array],DELIMITER,"-append");
dlmwrite("resultados_octave.txt",ro_array',DELIMITER,"-append");
endif

return;
