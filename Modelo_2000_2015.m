clear all
close all
clc
tStart = tic;
%% Modificado Por: Duilio Fonseca, Universidad Austral de Chile              
% Bajo el marco del proyecto de Balance Hídrico DGA.                        
% Modificaciones:                                                           
% El código requiere de un archivo shape (inventario glaciar) modificado,
% que considera en sus columnas la la suma, media y desviación              
% estándar de la diferencia de altura en tiempo para cada glaciar (dh/dt)   
% contacto: duilio.fonseca@uach.cl ; mschaefer@uach.cl     

%% Ajustes para corrida
%%% GLACIARES PARA PRUEBA
%%% Todos           (0) 
%%% Parinacota      (1)    
%%% Tronquitos      (2)    
%%% San Francisco   (3)    
%%% Bello           (4)    
%%% Mocho           (5)    
%%% Exploradores    (6)    
%%% Tyndall         (7)    
%%% Echaurren       (8)    
Glaciares_prueba=table();
Glaciares_prueba.codigo=[101010001;103420003;105702004;105703058;110301001;111421023;112288058;111440119];
Glaciares_prueba.nombre={'PARINACOTA';'TRONQUITOS';'SANFRANCISCO';'BELLO';'MOCHO';'EXPLORADORES';'TYNDALL';'OUTLIER'};

compu='Geofis';   %Computador de modelamiento
model='/V25_TcDEM_bandas';  %Nombre de corrida
%prefijo='Tc0_Tsr_hugonnet'; %Prefijo añadido al nombre de cada archivo
glaciar_prueba=0; %Glaciar o glaciares por modelar.
fr=0.2; % Factor multiplicativo para gl. de roca
dato_val=19; % 19 mean hugonnet 24 mean bandas hugonnet;
%% Forzar o automatizar zona glaciológica
%%% Condicionará método de modelación (Radiation-Temperature Index y Degree Day)
%%% Zona_var=0 %Forzar zona (ZN,ZC,ZS,ZA),
%%% Zona_var=1 %Selección de zona automática
zona_var=1;         % Método de selección de Zona
zona_forzada='ZS';  % Zona por forzar.

%% Set de iteración & Validación
%%% Para utilziar Braun et al. :braun
%%% Para utilizar Hugonnet et al. :hugonnet
%%% Preselección: método de validación para busqueda de mejor resultado para
%%% set de datos de Ts
%%% Selección: Elección de el mejor resultado segun set de datos de DDI/DDS o
%%% C1 y C0
pre_metodo='hugonnet';
metodo='hugonnet';
%pre_metodo='hugonnet';
%metodo='hugonnet';

%%% Rango de temperaturas para umbral de lluvia/nieve
Ts=(0:1:5); % [C] Temepratura para umbral de lluvia/nieve
%Ts=0;

% Switch para uso de lapse rate:
%%% (0)  Usar Tc_input "manual"
%%% (1)  Usar Tc de DEM;
%%% (2)  usar Lapserate de temperatura para corrección
lr_switch=1;    %SWITCH
Tc_input=0;     %T DE CORRECIÓN [C], manual ej: -2

% Periodo de modelación [Disponibilida de forzantes]
i_yr=2000;  % Año inicial de modelación
i_mth=3;    % Mes inicial de modelación
i_day=16;   % Dia inicial de modelación
f_yr=2019;  % Año final de modelación
f_mth=12;   % Mes final de modelación
f_day=31;   % Día final de modelación

%% 1. Direccionamiento y Periodo de modelamiento
%% SET GEOFIS
if compu=="Geofis" %Geofis
    
disp('Sección:1')
% 1.1 Asignación de directorio para resultados y nombre a modelación en particular
main='/home/duiliof/Desktop/BH/resultados'; % Establece dirección principal de resultados

% 1.2. Dirección de información geodésica para cada glaciar <REASIGNAR DIRECTORIO>.
link_gl = '/home/duiliof/Desktop/BH/input/info_plots';
dir_gl=dir(link_gl);
dir_gl=dir_gl(3:end);

% 1.3. Códigos de cada glaciar <NO MODIFICAR>.
codigos = dir(link_gl);
codigos = codigos(3:end);

% 1.4. Se cargan vectores de latitud y longitud para referenciar
% coordenadas <REASIGNAR DIRECTORIO>.
load /home/duiliof/Desktop/BH/input/geodesico_referencial/lat
load /home/duiliof/Desktop/BH/input/geodesico_referencial/lon

% 1.5. se carga la máscara de Chile para seleccionar solamente los pixeles
% dentro de la zona de estudio <REASIGNAR DIRECTORIO>.
load /home/duiliof/Desktop/BH/input/geodesico_referencial/Chile_mask.mat

% 1.6. link donde se encuentra la información meteorológica <REASIGNAR DIRECTORIO>.
link_meteo=strcat('/home/duiliof/Desktop/BH/input/Forcings_Glaciers_v25');

% 1.7. se carga el DEM para cálculo de diferencia de cota <REASIGNAR DIRECTORIO>.
[DEM header] = geotiffread('/home/duiliof/Desktop/BH/input/geodesico_referencial/Chile_DEM_005x005grad.tif');
DEM(DEM < 0) = NaN;

%%  SET Marius

elseif compu=="marius"
disp('Sección:1')
% 1.1 Asignación de directorio para resultados y nombre a modelación en particular
main='/home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/ResultadosModelacion'; % Establece dirección principal de resultados

% 1.2. Dirección de información geodésica para cada glaciar <REASIGNAR DIRECTORIO>.
link_gl = '/home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/Input/info_pruebas_HBH_tcmean';
dir_gl=dir(link_gl);
dir_gl=dir_gl(3:end);

% 1.3. Códigos de cada glaciar <NO MODIFICAR>.
codigos = dir(link_gl);
codigos = codigos(3:end);

% 1.4. Se cargan vectores de latitud y longitud para referenciar
% coordenadas <REASIGNAR DIRECTORIO>.
load /home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/Input/geodesico_referencial/lat
load /home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/Input/geodesico_referencial/lon

% 1.5. se carga la máscara de Chile para seleccionar solamente los pixeles
% dentro de la zona de estudio <REASIGNAR DIRECTORIO>.
load /home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/Input/shapes/Chile_mask.mat

% 1.6. link donde se encuentra la información meteorológica <REASIGNAR DIRECTORIO>.
link_meteo=strcat('/home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/Input/forzantes_25_setprueba');

% 1.7. se carga el DEM para cálculo de diferencia de cota <REASIGNAR DIRECTORIO>.
[DEM header] = geotiffread('/home/marius/Dropbox/Work/Varios/BalanceHidrico/BH5/Input/dem/Chile_DEM_005x005grad.tif');
DEM(DEM < 0) = NaN;

%% SET OMEN
elseif compu=="OMEN"
disp('Sección:1')

% 1.1 Asignación de directorio para resultados y nombre a modelación en particular
main='/home/duilio/Desktop/balance_hidrico-main/Resultados'; % Establece dirección principal de resultados

% 1.2. Dirección de información geodésica para cada glaciar <REASIGNAR DIRECTORIO>.
link_gl = '/home/duilio/Desktop/balance_hidrico-main/inputs/info_HBH';
dir_gl=dir(link_gl);
dir_gl=dir_gl(3:end);

% 1.3. Códigos de cada glaciar <NO MODIFICAR>.
codigos = dir(link_gl);
codigos = codigos(3:end);

% 1.4. Se cargan vectores de latitud y longitud para referenciar
% coordenadas <REASIGNAR DIRECTORIO>.
load /home/duilio/Desktop/balance_hidrico-main/inputs/geodesico_referencial/lat
load /home/duilio/Desktop/balance_hidrico-main/inputs/geodesico_referencial/lon

% 1.5. se carga la máscara de Chile para seleccionar solamente los pixeles
% dentro de la zona de estudio <REASIGNAR DIRECTORIO>.
load /home/duilio/Desktop/balance_hidrico-main/inputs/Chile_mask.mat

% 1.6. link donde se encuentra la información meteorológica <REASIGNAR DIRECTORIO>.
link_meteo=strcat('/home/duilio/Desktop/balance_hidrico-main/inputs/forzantes_25_setprueba');

% 1.7. se carga el DEM para cálculo de diferencia de cota <REASIGNAR DIRECTORIO>.
[DEM header] = geotiffread('/home/duilio/Desktop/balance_hidrico-main/inputs/Chile_DEM_005x005grad.tif');
DEM(DEM < 0) = NaN;

end

% 1.8. Creación de directorios para almacenamiento de resultados
% Establece dirección principal caudales por glaciar <NO MODIFICAR>.
mkdir(strcat(main,model))
output_link=strcat(main,model);

% 1.9. Creación de directorios <NO MODIFICAR>.
mkdir(strcat(output_link,'/Figuras_DD'));
mkdir(strcat(output_link,'/Qg_y_params_DD'));
mkdir(strcat(output_link,'/Resultados_por_glaciar'));
mkdir(strcat(output_link,'/DGA'));
mkdir(strcat(output_link,'/Figuras_ciclos'));
mkdir(strcat(output_link,'/input_plots'));
% 1.9. Directorio para guardar archivos de salida <NO MODIFICAR>
link_save2=strcat(output_link,'/Qg_y_params_DD');
link_save3=strcat(output_link,'/Resultados_por_glaciar');
link_rpublico=strcat(output_link,'/DGA');
link_ciclos=strcat(output_link,'/Figuras_ciclos');
link_inputplot=strcat(output_link,'/input_plots');
% 1.10 Directorio para las figuras de salida <NO MODIFICAR>
link_fig=strcat(output_link,'/Figuras_DD');
link_fig2=strcat(output_link,'/Figuras_forced');
mkdir(strcat(output_link,'/fichas'))

% 1.11 Ajustes de vectores de tiempo tiempo: datenum(año, mes, dia, hora, minuto, segundo)

% [periodo 1981-2015] tiempo_datenum = datenum(1981,1,1,1,0,0) : 1/8 : datenum(2020,4,30,0,0,0);
% Periodo de modelación validable 2000 a 2019.
tiempo_datenum = datenum(i_yr,i_mth,i_day,0,0,0) : 1/8 : datenum(f_yr,f_mth,f_day,0,0,0);
tiempo_completo_diario = datenum(1981,1,1):1:datenum(2020,4,30);
tiempo_completo_diario=tiempo_completo_diario';
tiempo = datevec(tiempo_datenum);
tiempo =[tiempo_datenum' tiempo];
clear tiempo_datenunm
tiempo_day = datenum(i_yr,i_mth,i_day) : 1 : datenum(f_yr,f_mth,f_day); %Periodo acotado para modelacion
tiempo_day = tiempo_day';
tiempo_day_vec=datevec(tiempo_day);
yr=tiempo_day_vec(:,1);
mth=tiempo_day_vec(:,2);
day=tiempo_day_vec(:,3);
filas_day = [];
clear tiempo

% indices para selección de datos en rango temporal validable.
index_2000=find(tiempo_completo_diario==datenum(2000,03,16));
index_2015=find(tiempo_completo_diario==datenum(2015,03,16));
index_2019=find(tiempo_completo_diario==datenum(2019,12,31));

% for i = 1:length(tiempo_day)
%     filas_day = [filas_day find((tiempo_day(i,1)) ==tiempo(:,2) & tiempo_day(i,2) == tiempo(:,3) & tiempo_day(i,3) == tiempo(:,4))];
% end

%% 2. Parámetros físicos para la modelación
disp('Sección 2: PARAMETROS')
% 2.1 Constante para área inicial

factor_Area = 0.5:0.1:1.5;
b=1.3598;             % Factor de ajuste Volume Area Scaling
c=0.027;              % Factor de ajuste Volume Area Scaling

% 2.2.  Constante de Stefann-Boltzmann
sigma = 5.67*10^(-8); % [ W / (m2 * K4)]

% 2.3  Emisividad de la nieve. Rango: 0.97 a 1.00
snow_emissivity = 0.99;

% 2.4. Emisividad del aire.
air_emissivity = 1;

% 2.5. Temperatura de la superficie del hielo y manto
T_snow = -1;

% 2.6. albedo máximo y mínimo que de la nieve y el glaciar
alb_max_snow = 0.90;
alb_min_snow = 0.60;
alb_ice = 0.40;       % Albedo del hielo
% 2.7. Temperatura de separación lluvia/nieve
To= 0; % [C] Temperatura para umbral de derretimiento

%% 2.7.1  Lapse Rate 
% 2.8 Switch para uso de lapse rate: (Muteado para pruebas)
% (0)  Usar Tc_input "manual"
% (1)  Usar Tc de DEM;
% (2)  usar Lapserate de temperatura para corrección
% lr_switch=1;    %SWITCH
% Tc_input=0;     %T DE CORRECIÓN

% 2.9 Parametrización para glaciares de zona Norte
%C1_v= 4:2:12;   % Rango para C1
%C0_v=-50:10:-10; % Rango para C0
C1_v= 2:2:12;   % Rango para C1
C0_v=-80:10:-10; % Rango para C0
C2_v=[0.6 0.8 1.0]; % Rango SW
% 2.10 Decaimiento del albedo
r = 0.10; % [1/d]
glaciares_no_modelados = []; %Creación de matriz vacia para modelados fallidos

% 2.11 Parametrización para balance de energía
k_o=0.4;                    %constante de von karmann
z= 2.0;                     %Altura de forzante
z_o= 0.0005;                % Parametro de rugosidad (extraido de Cuffey and paterson 2010)
C=((k_o^2)/log(z/z_o)^2);   % C*  Coeficiente de transferencia efectiva
t_s=0;                      %Temperatura superficial °C
Avps=611;                   %¨Presion de vapor en la superificie [Pascal]

%% 3.Preparar DDI, DDS, C1 y C2: estos parámetros son las únicas variables que jugaran en el seteo definitivo.
DD_I_i=1:1:14;       %Rango de DDF de Hielo [mm/C día] (1:11 primera corrida)
DD_S_i=0.1:0.25:0.9;   %Rango porcentual de  DDF de nieve [mm/C día] (1:11 primera corrida)

% 3.1.1.  Creación de vectores con Posibles pares de DDI y DDS
for hh=1:length(DD_I_i)
    for rr=1:length(DD_S_i)
        DD_S_ii(rr,hh)=DD_I_i(hh).*DD_S_i(rr);           %porcentaje
        DD_I_ii(:,hh)=DD_I_i(hh)*ones(length(DD_S_i),1);
    end
end
clear DD_I_i DD_S_i
DD_I=reshape(DD_I_ii,[],1);
DD_S=reshape(DD_S_ii,[],1);

pairs=[DD_I(:) DD_S(:)];        %Matriz de posibles pares de DDI y DDS
clear DD_I DD_S

% 3.4. Calculo de pares de C1 y C2 para glaciares de Zona Norte y Zona
% Central <NO MODIFICAR>
[C1,C0] = meshgrid(C1_v,C0_v);
prepairsz=cat(2,C1',C0');
pairs_zn_pre=reshape(prepairsz,[],2);
clear C1_v C0_v C1 C0 prepairsz
for cc=1:length(C2_v)
C2(:,cc)=ones(length(pairs_zn_pre),1)*C2_v(cc);
end
pairs_zn=zeros(length(pairs_zn_pre)*3,3);
for i=1:size(C2,2)
pairs_zn((length(C2)*(i-1)+1):length(C2)*(i),:)=[pairs_zn_pre C2(:,i)];
end
%% 4. Punto de inicio del modelo
disp('Sección 3: Carga de glaciares')

% 4.1 punto de partida de iteración (equivalente a glaciar puntual desde el cual se inicia la iteración)
 l_inicio=1;

% 4.2 Punto Final de iteración, equivalente al último glaciar por iterar
l_final = length(dir_gl);


% 4.3. Loop principal: la longitud de éste corresponde a la cantidad de archivos
% a leer.
for i = l_inicio:l_final
 %for i=42 %limitado a glaciares de Zona Norte

    % 4.4 Prueba para un glaciar puntual (comentar bloque 4.3 y descomentar 4.4).
    % for i=16
    %   codigo=codigos{i,1};
    %    codigo=codigo(3:end);

    % 4.5. se carga la data geodésica de cada glaciar, se verifica que cuente con dh/dt, de no
    % ser así pasa al segundo loop de iteración.
   
  %%% [ESTRUCTURA DEL INPUT] Codigo(1) Latitud(2) Longitud(3) AÑO(4) MES(5) DIA(6) Area(7) Areas2(8) Volumen(9)
  %%% Pendiente(10) Orientacion(11) altura_promedio(12) tipo(13) Braun_dhdt(SUM)(14)
  %%% Braun_dhdt(MEAN)(15) Braun_dhdt(SD)(16) Braun_dhdt(COUNT)(17) Hugonnet_dhdt(SUM)(18)
  %%% Hugonnet_dhdt(MEAN)(19) Hugonnet_dhdt(SD)(20) Hugonnet_dhdt(COUNT) (21)
  %%% Calving_factor(gt/y) (22)
  
   %load(strcat(link_gl,'/GLA_CL',num2str(Glaciares_prueba.codigo(i)),'.mat'))
   load(strcat(link_gl,'/',codigos(i).name))
   %load(strcat(link_gl,'/',codigos(exp_no_mod(i)).name))
    try
    Tc_dem=info(23);
    catch 
    Tc_dem=0;
    end 
    
    % 4.5.1 En la celda 14 del archivo info, se encuentra el valor extraido del
    %producto de Braun et al. (2018)
    %¿Cuenta el glaciar con información para validación post modelado?
    if info(1,14)~=0 || info(1,18)~=0 
        ok=true;	%si
    else
        ok=false;	%no
    end

    %% 5. Cálculo de datos geodésicos y georeferenciales
    disp('Sección 4: Calculo de datos geodésicos')
    
    % 5.1 se encuentra la posición de la latitud y longitud más cercana del glaciar respecto a la resolución de grilla de NetCDF (Cr2Met)
    fila_2000 = find(abs(lat-info(2)) == min(abs(lat-info(2))));
    col = find(abs(lon-info(3)) == min(abs(lon-info(3))));
    lat1=lat(fila_2000);
    lon1=lon(col);

    % 5.1.1 Clasificación de Macrozona [Nuevo]
    %Zona Norte 18°-32°
    if zona_var==1 %Añadido solo para pruebas [quitar para version final]
        if lat1>-32 && lat1<-16
            zona='ZN';
            %Zona Centro 32°- 36°
        elseif lat1>-36 && lat1<-32
            zona='ZC';
            %Zona Sur 36°-46°56°
        elseif lat1>-46 && lat1<-36
            zona='ZS';
            %Zona Austral 46°-
        elseif lat1>-56 && lat1<-46
            zona='ZA';
        end
    elseif zona_var==0 %Añadido solo para pruebas [quitar para version final]
        zona=zona_forzada; %Añadido solo para pruebas [quitar para version final]
    end %Añadido solo para pruebas [quitar para version final]

    % 5.2. Calculo de diferencia de cota entre glaciar y elev. media del glaciar (no tiene
    % relevancia en el código.
    diff_cota = DEM(fila_2000, col) - info(1,12);
    diff_cota=double(diff_cota);

    % 5.3 Variación de volumen 2000-2011/2015 Braun et al. Se multiplica la suma de
    % dh/dt por el área de pixel (900m^2) y se convierte a km^3
    
    % Calculo del área del pixel mediante la formula de Haversine
    res_x=0.000277777800000000429; %resolusión pixel grados producto braun et al.
    res_y=0.000277777800009570031; %resolusión pixel grados producto braun et al.
    [area_px_braun]=areapx(lat1,lon1,res_x,res_y); % Area por pixel [m2]
    clear res_x res_y
    
    res_x=0.0009008701491566568735; %resolusión pixel grados producto braun et al.
    res_y=0.0009008701491602962149; %resolusión pixel grados producto braun et al.
    [area_px_ht]=areapx(lat1,lon1,res_x,res_y); % Area por pixel [m2]
    clear res_x res_y
    
    %% 7. Cargando datos meteorológicos
    disp('Sección 5: INICIO DE LOOP')
    % 7.0. Corrección mediante gradiente de temperatura -6.49 K/km ISA (international standar atmosphere) 0-11000 msnm
    if lr_switch==2
        lr=6.49/1000;
        Tc=(diff_cota*lr);
    elseif lr_switch==1
        Tc=Tc_dem;
    elseif lr_switch==0
        Tc=Tc_input;
    end

    % 7.1 Verificar si existe data meterológica en relación al glaciar modelado.

    if exist(strcat(link_meteo,'/data_',num2str(lat(fila_2000)),'_',num2str(lon(col))),'file') > 0 || exist(strcat(link_meteo,'/data_CL',num2str(info(1)),'.txt'),'file')>0 %Dato
        try %574

            % 7.2 Si el glaciar tiene un área mayor a 25 km^2 se carga temperatura de promedio de varios puntos
            if info(:,7)/10^6>=25
                file_temp_poligonal=readmatrix(strcat(link_meteo,'/data_CL',num2str(info(1)),'.txt'));
                temperatura=file_temp_poligonal(:,7);
                pr=file_temp_poligonal(:,2);
                ws=file_temp_poligonal(:,3);
                pa=file_temp_poligonal(:,4);
                sw=file_temp_poligonal(:,5);
                lw=file_temp_poligonal(:,6);
                t2m=file_temp_poligonal(:,7)+Tc;
                rh=file_temp_poligonal(:,8);


                % 7.3. Revisión de input de temperatura
                if isnan(temperatura(1))==false
                    % 7.4. Si el archivo no contiene información se carga el dato de temperatura a partide un solo punto.
                elseif exist('temperatura')==false
                    file_temp_puntual=dlmread(strcat(link_meteo,'/data_',num2str(lat(fila_2000)),'_',num2str(lon(col))),'\t');
                    pr=file_temp_puntual(:,1);
                    ws=file_temp_puntual(:,2);
                    pa=file_temp_puntual(:,3);
                    sw=file_temp_puntual(:,4);
                    lw=file_temp_puntual(:,5);
                    t2m=file_temp_puntual(:,6)+Tc;
                    rh=file_temp_puntual(:,7);

                    % 7.5. Corrección mediante gradiente de temperatura -6.49 K/km ISA (international standar atmosphere) 0-11000 msnm

                end % En de if (231)
            else
                file_temp_puntual=load(strcat(link_meteo,'/data_',num2str(lat(fila_2000)),'_',num2str(lon(col))),'\t');
                 pr=file_temp_puntual(:,1);
                ws=file_temp_puntual(:,2);
                pa=file_temp_puntual(:,3);
                sw=file_temp_puntual(:,4);
                lw=file_temp_puntual(:,5);
                t2m=file_temp_puntual(:,6)+Tc;
                rh=file_temp_puntual(:,7);
                pr=sum_timeseries(pr,8);
                ws=average_timeseries(ws,8);
                pa=average_timeseries(pa,8);
                sw=average_timeseries(sw,8);
                lw=average_timeseries(lw,8);
                t2m=average_timeseries(t2m,8);
                rh=average_timeseries(rh,8);

            end

            %% 8. Carga de áreas, dv/dt y da/dt
            disp('Sección 6: Datos geodésicos')
            if ok==true %Condicionantes: Si el glaciar cuenta o no con dato de dh/dt
                
                % Cálculo de área-volumen, dV/dt, Variación de área estimada mediante el producto de Braun et al.
                dvdt_braun_2000_2014=info(1,14)*area_px_braun/10^9; % [km^3]
                dvdt_hugonnet_2000_2019=info(1,18)*area_px_ht/10^9; % [km^3]  
                
                % Estimación de volumen observado en inventario. ESTO VA A
                % GENERAR UN PROBLEMA CON EL NUEVO INVENTARIO.
                v_obs=c*(info(:,7)/10^6)^b; %[km^3]
                a_obs=(info(:,7)/10^6);
                a_obs_ac(:,i)=a_obs;
                
                % Se diferencia si hay un solo año con información o si
                % Calcular la diferencia entre la observación y la
                % estimación
                t_obs = datenum(info(6),03,16);
                
                % Volumenes y áreas estimadas mediante Braun et al.                
                V_2014_est_braun=v_obs+dvdt_braun_2000_2014*(2014-info(:,6));  % Volumen estimado para 2015 a partir de producto de Braun et al.
                A_2014_est_braun=(V_2014_est_braun/c)^(1/b); 
                        
                % Volumenes y áreas estimadas mediante Hugonnet.
                V_2019_est_ht=v_obs+dvdt_hugonnet_2000_2019*(2019-info(:,6));  % Volumen estimado para 2015 a partir de producto de Braun et al.
                A_2019_est_ht=(V_2019_est_ht/c)^(1/b);
                                    
                
                % Si el área del inventario es muy antigua se asume como 20
                if t_obs<datenum(2000,03,16,00,00,00)
                t_obs=datenum(2000,03,16,00,00,00);
                else
                    %nada
                end
                
                %Vector de tiempo de la modelación
                t_mod=tiempo_day(find(tiempo_day==t_obs)):1:tiempo_day(end); %Obtención de periodo validable
                t_mod_vec=datevec(t_mod);
                t_mod_yr=t_mod_vec(:,1);
                t_mod_mth=t_mod_vec(:,2);
                t_mod_day=t_mod_vec(:,3);
                t_mod=t_mod';
                yr_mod=tiempo_day_vec(find(tiempo_day==t_obs),1);
                time_plot=datetime(t_mod,'ConvertFrom','datenum');
                %Creación de indices para gŕaficos
                index_mod_i=find(tiempo_completo_diario==t_mod(1));
                index_2014_plot=find(t_mod==datenum(2014,03,16));
                index_2019_plot=find(t_mod==datenum(2019,03,16));
                %fila_obs = find(tiempo_day == t_obs);
                fila_obs=1;
                
                
                                % Selección de zona
                if zona=="ZN" || zona=="ZC"
                    disp('Sección 7(a): Modelación RTI ')
                    % Iteración para la modelación con distintos DDI y DDS.
                    for aa=1:length(pairs_zn)

                        vapourpressuresaturation=(6.116441*10.^((7.591386.*t2m)./(t2m+240.7263))).*100;
                        vapourpressure=(rh.*vapourpressuresaturation);

                        % Input para modelación de derretimiento
                        % [Version completa] x_diario=[tiempo_day pr t2m rh pa ws vapourpressure sw lw];
                        x_diario=[t_mod pr(index_mod_i:index_2019) t2m(index_mod_i:index_2019) rh(index_mod_i:index_2019) pa(index_mod_i:index_2019) ws(index_mod_i:index_2019) vapourpressure(index_mod_i:index_2019) sw(index_mod_i:index_2019) lw(index_mod_i:index_2019)]; % Versión acotada
                   
                    for qq=1:length(Ts)
                        if info(1,13) == 1 % verificar si el glaciar es rocoso
                            [parametros, flujos] = balance_ZN_tiempo_acotado(x_diario, To, Ts(qq), a_obs, c, b, pairs_zn(aa,:), alb_ice, fr, C, Avps,t_s,yr_mod ,info(22));
                        else
                            [parametros, flujos] = balance_ZN_tiempo_acotado(x_diario, To, Ts(qq), a_obs, c, b, pairs_zn(aa,:), alb_ice, 1, C, Avps,t_s,yr_mod ,info(22));
                        end
                    flujos_acum_pre(:,:,qq)=flujos;         %Resultados previos para cada Ts posible
                    parametros_acum_pre(:,:,qq)=parametros; %Resultados previos para cada Ts posible
                    end
                    
                    %Selección de mejor resultado de cada set de Ts
                   [bn_nmin_dhdt_pre,bn_dhdt_modelado_pre,bn_dhdt_diferencia_pre,bn_dhdt_medido_pre,bn_dhdt_diferencia_r_pre]=validacion_dhdt_braun(t_mod,info(15),flujos_acum_pre,yr_mod);
                   [ht_nmin_dhdt_pre,ht_dhdt_modelado_pre,ht_dhdt_diferencia_pre,ht_dhdt_medido_pre,ht_dhdt_diferencia_r_pre]=validacion_dhdt_hugonnet(t_mod,info(dato_val),flujos_acum_pre,yr_mod);
                   
                   
                   if pre_metodo=="braun"
                        pre_seleccion=bn_nmin_dhdt_pre;  
                         if isempty(pre_seleccion)==1
                                break
                            else
                           %nada        
                         end
                        
                   elseif pre_metodo=="hugonnet"
                    pre_seleccion=ht_nmin_dhdt_pre;  
                        if isempty(pre_seleccion)==1
                                break
                            else
                           %nada        
                         end
                   end
                   
                                    
                   % Método de preseleccion segun distintos MB geodésicos
                   % Si existe más de un mínimo tomar el primer valor
                   if length(pre_seleccion)>1
                    pre_seleccion=pre_seleccion(1);
                   else
                        %nada
                   end

                   %pre_seleccion=bn_nmin_dhdt_pre;
                   % Almacenamiento de resultados posibles 
                   flujos_acum(:,:,aa)=flujos_acum_pre(:,:,pre_seleccion);
                   parametros_acum(:,:,aa)=parametros_acum_pre(:,:,pre_seleccion);
                    
                    % Reset de parametros
                    clear bn_nmin_dhdt_pre bn_dhdt_modelado_pre bn_dhdt_diferencia_pre bn_dhdt_medido_pre flujos_acum_pre
                    clear ht_nmin_dhdt_pre ht_dhdt_modelado_pre ht_dhdt_diferencia_pre ht_dhdt_medido_pre parametros_acum_pre
                    end %End for de pirs_zn


                elseif zona=="ZS" || zona=="ZA"
                    disp('Sección 7(b): Modelación DDI ')
                    %% Solo para pruebas manuales
                    % pairs=[3.5,0.35];
                    % for ee=1:1
                    for ee=1:length(pairs)

                        % Input para modelación de derretimiento
                        x_diario=[t_mod pr(index_mod_i:index_2019) t2m(index_mod_i:index_2019)];

                   for qq=1:length(Ts)
                        if info(1,13) == 1 % verificar si el glaciar es rocoso
                            [parametros, flujos] = balance_ZA_tiempo_acotado(x_diario, To, Ts(qq),  a_obs, c, b, pairs(ee,:), alb_ice, fr, yr_mod, info(22));
                        else
                            [parametros, flujos] = balance_ZA_tiempo_acotado(x_diario, To, Ts(qq), a_obs, c, b, pairs(ee,:), alb_ice, 1, yr_mod, info(22));
                        end
                    flujos_acum_pre(:,:,qq)=flujos;         %Resultados previos para cada Ts posible
                    parametros_acum_pre(:,:,qq)=parametros; %Resultados previos para cada Ts posible
                   end
                   
                   %Selección de mejor resultado de cada set de Ts
                   [bn_nmin_dhdt_pre,bn_dhdt_modelado_pre,bn_dhdt_diferencia_pre,bn_dhdt_medido_pre,bn_dhdt_diferencia_r_pre]=validacion_dhdt_braun(t_mod,info(15),flujos_acum_pre,yr_mod);
                   [ht_nmin_dhdt_pre,ht_dhdt_modelado_pre,ht_dhdt_diferencia_pre,ht_dhdt_medido_pre,ht_dhdt_diferencia_r_pre]=validacion_dhdt_hugonnet(t_mod,info(dato_val),flujos_acum_pre,yr_mod);
                   
                   % Método de preseleccion segun distintos MB geodésicos
                   if pre_metodo=="braun"
                    pre_seleccion=bn_nmin_dhdt_pre;
                          if isempty(pre_seleccion)==1
                                break
                            else
                           %nada        
                         end
                   elseif pre_metodo=="hugonnet"
                    pre_seleccion=ht_nmin_dhdt_pre;  
                         if isempty(pre_seleccion)==1
                                break
                            else
                           %nada        
                         end
                   end 
                   
                   % Si existe más de un mínimo tomar el primer valor
                   if length(pre_seleccion)>1
                    pre_seleccion=pre_seleccion(1);
                   else
                        %nada
                   end

                    
                   % Almacenamiento de resultados posibles                       
                   flujos_acum(:,:,ee)=flujos_acum_pre(:,:,pre_seleccion);
                   parametros_acum(:,:,ee)=parametros_acum_pre(:,:,pre_seleccion);
                   
                   % Reset de parametros
                   clear bn_nmin_dhdt_pre bn_dhdt_modelado_pre bn_dhdt_diferencia_pre bn_dhdt_diferencia_pre_r bn_dhdt_medido_pre flujos_acum_pre
                   clear ht_nmin_dhdt_pre ht_dhdt_modelado_pre ht_dhdt_diferencia_pre ht_dhdt_diferencia_pre_r ht_dhdt_medido_pre parametros_acum_pre

                    end %end del for de pairs
                end %end para condicional de zonas
                             
               disp('Sección 8: VALIDACION')           
               
               %% Validación experimental con MAPE [RF]
               % [MAPE_fit_best,MAPE_raw_best,MAPE_fit_index,MAPE_raw_index]=validacion_MAPE(v_obs,var_volumen_2000_2015,flujos_acum);
               % [RMSE_fit_best,RMSE_raw_best,RMSE_fit_index,RMSE_raw_index]=validacion_RMSE(v_obs,var_volumen_2000_2015,flujos_acum);
               
               [bn_nmin_dhdt,bn_dhdt_modelado,bn_dhdt_diferencia,bn_dhdt_medido,bn_dhdt_diferencia_r]=validacion_dhdt_braun(t_mod,info(15),flujos_acum,yr_mod);
               % Validación mediante balance de masa contrastado con braun.
               [ht_nmin_dhdt,ht_dhdt_modelado,ht_dhdt_diferencia,ht_dhdt_medido,ht_dhdt_diferencia_r]=validacion_dhdt_hugonnet(t_mod,info(dato_val),flujos_acum,yr_mod);
                
               %% Selección de mejor simulación
               % Otros sistemas de seleccion
                %puntaje_menor=MAPE_raw_index;
                %puntaje_menor=RMSE_raw_index;
                %puntaje_menor=find(abs(puntaje)==min(abs(puntaje)));
                
                %Metodos de seleccion [Dejado en set inicial para pruebas]
                 if metodo=="braun"
                   seleccion=bn_nmin_dhdt; %Selección validad mediante Braun et al.
                 elseif metodo=="hugonnet"
                   seleccion=ht_nmin_dhdt; %Selección validad mediante Braun et al.
                 end 
                
                    
                %seleccion=bn_nmin_dhdt; %Selección validad mediante Hugonett et al.
                
                %9.8.7. Si todas las simulaciones obtienen un puntaje equivalente el modelo toma el primer resultado. (review)
                
                if length(seleccion)>1
                    seleccion=seleccion(1);
                else
                    %nada
                end

                % 9.8.8. Seleccion de mejor simulacion a partir de dataset.
                parametros_gl=parametros_acum(:,:,seleccion);
                flujos_gl=flujos_acum(:,:,seleccion);
                clear parametros_acum flujos_acum
                
            end  %If ok==true

            % 9.8.9. Si existen resultados, se ordena continuar con la sección de graficado y creación de archivos.
            if exist('flujos_gl')==1
                % Continua al gráficado
                continuar=true;

                % 9.8.10 Hay resultados, pero no superan el primer threshold por consiguiente pasa
                % a ser modelado en el modulo original
            elseif exist('flujos')==1 & exist('flujos_gl')==0
                continuar=false;
                tiene_dhdt=true;
                clear flujos parametros

                % 9.8.11 No existen resultados se modela en modulo original
            elseif exist('flujos_gl')==0
                continuar=false;
                clear flujos parametros
            end
            
        
            if continuar==true
                %% Calculo de Balance de masa, ciclos para forzantes y creación de tablas
                 % 10.5 Creación de tablas Para DGA y publico general
                if zona=="ZN" || zona=="ZC"
                    Cabeceras_diarias = {'Año';'Mes';'Día';'Temperatura [°C]';'Precipitación [mm]';'Precipitación sólida[mm]';'Área[km2]';'Volumen[km3]';'SWE [mm]';'Derretimiento nieve[mm]';'Derretimiento hielo [mm]';'Caudal Nival[m3/s]';'Caudal Hielo[m3/s]';'Caudal Total[m3/s]'};
                    tabla_diaria=table(t_mod_yr,t_mod_mth,t_mod_day,t2m(index_mod_i:index_2019),pr(index_mod_i:index_2019),flujos_gl(:,13),flujos_gl(:,2),flujos_gl(:,3),flujos_gl(:,8),flujos_gl(:,9),flujos_gl(:,7),flujos(:,14),flujos(:,12),flujos(:,15),'VariableNames',Cabeceras_diarias);
                   % Cabeceras_info={'Codigo Glaciar';'Latitud [°]';'Longitud [°]';'Altitud [msnm]';'Área inventario [km2]';'Zona glaciológica';'Rocoso [boolean]';'dh/dt[m/año]';'Suma de errores porc. [%]';'C1';'C0'};
                   % tabla_informacion_glaciar=table(info(1),info(2),info(3),info(1,12),(info(:,7)/10^6),string(zona),info(1,13),info(1,14),puntaje(puntaje_menor),parametros_gl(1,7),parametros_gl(1,8),'VariableNames',Cabeceras_info);
                elseif zona=="ZA" || zona=="ZS"
                    Cabeceras_diarias = {'Año';'Mes';'Día';'Temperatura [°C]';'Precipitación [mm]';'Precipitación sólida[mm]';'Área[km2]';'Volumen[km3]';'SWE [mm]';'Derretimiento nieve[mm]';'Derretimiento hielo [mm]';'Caudal Nival[m3/s]';'Caudal Hielo[m3/s]';'Caudal Total[m3/s]'};
                    tabla_diaria=table(t_mod_yr,t_mod_mth,t_mod_day,t2m(index_mod_i:index_2019),pr(index_mod_i:index_2019),flujos_gl(:,13),flujos_gl(:,2),flujos_gl(:,3),flujos_gl(:,8),flujos_gl(:,9),flujos_gl(:,7),flujos(:,14),flujos(:,12),flujos(:,15),'VariableNames',Cabeceras_diarias);
                   % Cabeceras_info={'Codigo Glaciar';'Latitud [°]';'Longitud [°]';'Altitud [msnm]';'Área inventario [km2]';'Zona glaciológica';'Rocoso [boolean]';'dh/dt[m/año]';'Suma de errores porc. [%]';'DDI';'DDS'};
                   % tabla_informacion_glaciar=table(info(1),info(2),info(3),info(1,12),(info(:,7)/10^6),string(zona),info(1,13),info(1,14),puntaje(puntaje_menor),parametros_gl(1,7),parametros_gl(1,8),'VariableNames',Cabeceras_info);
                end

                % 10.6 en caso de que la validación falle por Sobrederretimiento (vol=0 y área=0). 
                if isempty(ht_dhdt_diferencia)==1
                    ht_dhdt_diferencia=NaN;
                    ht_dhdt_diferencia_r=NaN;
                    ht_dhdt_modelado=NaN;
                    ht_nmin_dhdt=NaN;
                else 
                    % Nada    
                end 
                
                if isempty(bn_dhdt_diferencia)==1
                    bn_dhdt_diferencia=NaN;
                    bn_dhdt_diferencia_r=NaN;
                    bn_dhdt_modelado=NaN;
                    bn_nmin_dhdt=NaN;
                else
                    % Nada
                end 
                
                if lat1>-36 %1 Norte
                    std_data_for_map(i,:)=[info(1) 1 lon1 lat1 info(6) a_obs v_obs yr_mod  parametros_gl(1,7) parametros_gl(1,8) bn_dhdt_medido bn_dhdt_modelado bn_dhdt_diferencia bn_dhdt_diferencia_r ht_dhdt_medido ht_dhdt_modelado ht_dhdt_diferencia ht_dhdt_diferencia_r Tc parametros_gl(1,1) parametros_gl(1,end)];
                else lat1<=-36; % 2 Sur
                    std_data_for_map(i,:)=[info(1) 2 lon1 lat1 info(6) a_obs v_obs yr_mod  parametros_gl(1,7) parametros_gl(1,8) bn_dhdt_medido bn_dhdt_modelado bn_dhdt_diferencia bn_dhdt_diferencia_r ht_dhdt_medido ht_dhdt_modelado ht_dhdt_diferencia ht_dhdt_diferencia_r Tc parametros_gl(1,1) nan];
                end
                % 10.7 Escrituta de archivos xlsx para DGA
                 % warning('off','MATLAB:xlswrite:AddSheet'); %optional
                 % writetable(tabla_diaria,strcat(link_rpublico,'/CL_',num2str(info(1,1)),'_',string(zona),'_',prefijo,'.xlsx'),'Sheet',1);
                  %writetable(tabla_informacion_glaciar,strcat(link_rpublico,'/CL_',num2str(info(1,1)),'_',string(zona),'_',prefijo,'.xlsx'),'Sheet',2);

                % Creación de tablas para análisis de ciclos.             
                tabla_diaria.("Caudal Hielo[m3/s]")(tabla_diaria.("Caudal Hielo[m3/s]")==0)=NaN;
                tabla_diaria.("Caudal Nival[m3/s]")(find(tabla_diaria.("Caudal Nival[m3/s]")==0))=NaN;
                tabla_diaria.("Caudal Total[m3/s]")(find(tabla_diaria.("Caudal Total[m3/s]")==0))=NaN;
                tabla_diaria.("Precipitación [mm]")(find(tabla_diaria.("Precipitación [mm]")==0))=NaN;
                tabla_diaria.("Derretimiento hielo [mm]")(find(tabla_diaria.("Derretimiento hielo [mm]")==0))=NaN;
                tabla_diaria.("Derretimiento nieve[mm]")(find(tabla_diaria.("Derretimiento nieve[mm]")==0))=NaN;
                tabla_diaria.("Precipitación sólida[mm]")(find(tabla_diaria.("Precipitación sólida[mm]")==0))=NaN;
                % Año hidrológico (15 de marzo)
                index_hy=find(tabla_diaria.Mes==3 & tabla_diaria.("Día")==16 );
                for ff=1:length(index_hy)-1
                    for pp=index_hy(ff):index_hy(ff+1)
                           yr_hidro(pp,1)=ff+yr_mod-1;
                    end
                end
                
                % Corte de tabla hasta fecha límite de año hidrológico;
                tabla_diaria(length(yr_hidro)+1:end,:)=[];
                tabla_diaria.("Año Hidrológico")=yr_hidro;
                promedio_mensual=varfun(@nanmean,tabla_diaria,'GroupingVariables',{'Año' 'Mes'},'InputVariables',{'Temperatura [°C]' 'Volumen[km3]' 'Área[km2]' 'SWE [mm]' 'Caudal Nival[m3/s]' 'Caudal Hielo[m3/s]' 'Caudal Total[m3/s]'});
                suma_mensual=varfun(@nansum,tabla_diaria,'GroupingVariables',{'Año' 'Mes'},'InputVariables',{'Precipitación [mm]' 'Precipitación sólida[mm]' 'Derretimiento nieve[mm]' 'Derretimiento hielo [mm]'});
                suma_anual=varfun(@nansum,tabla_diaria,'GroupingVariables',{'Año Hidrológico'},'InputVariables',{'Precipitación [mm]' 'Precipitación sólida[mm]' 'Derretimiento nieve[mm]' 'Derretimiento hielo [mm]' 'SWE [mm]'});
                clear yr_hidro index_hy
                % Creación de tablas para MB
                MB=table();
                MB.yr=t_mod_yr;
                MB.mth=t_mod_mth;
                MB.day=t_mod_day;
                MB.area=flujos_gl(:,2);
                MB.volumen=flujos_gl(:,3);
                MB.dh=(flujos_gl(:,3)./flujos_gl(:,2))*1000; % [m.w.eq]
                MB_year=varfun(@diff,MB,'GroupingVariables',{'yr'});
                MB_year(find(MB_year.diff_dh==0),:)=[];
                               
                % Periodización de precipitaciones en rangos temporales
                % para ciclo anual
                idx(:,1)=(suma_mensual{:,1}>=yr_mod) & (suma_mensual{:,1}<2019);
                   
                %Extracción de datos
                pr_1=suma_mensual{idx(:,1),[1 2 4]};
                
                %Creacion de tablas vacias
                pr_idx1=table();
               
                %Asignacion de variables
                pr_idx1.('Año')=pr_1(:,1);
                pr_idx1.('Mes')=pr_1(:,2);
                pr_idx1.('Precipitaciones [mm]')=pr_1(:,3);
                clear pr_1

                %Calculo de ciclo anual para cada periodo
                pr_idx1=varfun(@mean,pr_idx1,'GroupingVariables',{'Mes'});
                
                datos_mensuales=join(promedio_mensual,suma_mensual);
                ciclo_anual=varfun(@nanmean,datos_mensuales,'GroupingVariables',{'Mes'});
                ciclo_anual_sd=varfun(@nanstd,datos_mensuales,'GroupingVariables',{'Mes'});
                ciclo_anual.Properties.VariableNames = {'Mes','GroupCount','Año','GroupCountsd','Temperatura [°C]','Volumen[km3]','Área[km2]','SWE [mm]','Caudal Nival[m3/s]','Caudal Hielo[m3/s]','Caudal Total[m3/s]','Precipitación sólida[mm]','Precipitación [mm]','Derretimiento nieve[mm]','Derretimiento hielo [mm]'};
                ciclo_anual_sd.Properties.VariableNames = {'Mes','GroupCount','Año','GroupCountsd','Temperatura [°C]','Volumen[km3]','Área[km2]','SWE [mm]','Caudal Nival[m3/s]','Caudal Hielo[m3/s]','Caudal Total[m3/s]','Precipitación sólida[mm]','Precipitación [mm]','Derretimiento nieve[mm]','Derretimiento hielo [mm]'};
                
                disp('Sección FIN: Escritura')      
                %Guardar archivo para input_plots
                save(strcat(link_inputplot,'/CL',num2str(info(1))),'time_plot','fila_obs','a_obs','v_obs','flujos_gl','info','index_2014_plot','index_2019_plot','index_mod_i','index_2019','A_2014_est_braun','A_2019_est_ht','V_2014_est_braun','V_2019_est_ht','t2m','sw','yr_mod','flujos_gl','parametros_gl','suma_anual','MB_year','bn_dhdt_modelado','bn_dhdt_medido','ht_dhdt_modelado','ciclo_anual','ciclo_anual_sd','pr_idx1') 
                
               
                % Guardado de archivos para equipo UCh y balance
                codigos_simulados(i,1)=info(1);
                %save(strcat(link_save2,'/Q_y_params_GL',num2str(info(1)),'_',string(zona),'_',prefijo), 'parametros_gl' , 'flujos_gl')
                save(strcat(link_save2,'/Q_y_params_GL',num2str(info(1))), 'parametros_gl' , 'flujos_gl')

                if exist('glaciares_incompetentes')==true
                    %  save('C:\Users\duili\Desktop\paper_nuevo\resultados\beagle\problemas','glaciares_incompetentes');
                else
                    %
                end

            end
            % 10.5. Se almacenan los glaciares que no pueden ser modelados
        catch
            glaciares_no_modelados = [glaciares_no_modelados, i];
            fprintf('Alerta de crash')
        end
        clc
        disp(strcat('Avance:', num2str(100*i/(l_final)),'%'))
    end
    % 10.6 Grabando datos para comparación de resultados
    
    clear puntaje puntaje_menor flujos_acum parametros_acum DD_I_m DD_S_m DD_I_r pr RMSE temperatura DD_S_r ok diferencia continuar  dif_vols_mod_ac x_diario dif_area_2015_mod_ac dif_area_obs_ac dif_area_obs_ac dif_braun_vol_1979_2000 dif_vols_acum diff_Area  diff_Area_obs DD_I DD_S flujos_gl parametros_gl clase;
    clear bn_nmin_dhdt bn_dhdt_modelado bn_dhdt_diferencia bn_dhdt_medido 
    clear ht_nmin_dhdt ht_dhdt_modelado ht_dhdt_diferencia ht_dhdt_medido 
    clear suma_mensual yr_mod i_yr i_mth i_day f_yr f_mth f_day idx info MB t_mod time_plot index_mod_i index_2014_plot index_2015_plot
    clear pr_idx1 ciclos suma_mensual promedio_mensual ciclo_anual ciclo_anual_sd tabla_diaria tabla_informacion_glaciar Cabeceras_diarias Cabeceras_info
    clear bn_nmin_dhdt_pre bn_dhdt_modelado_pre bn_dhdt_diferencia_pre bn_dhdt_medido_pre flujos_acum_pre
    clear ht_nmin_dhdt_pre ht_dhdt_modelado_pre ht_dhdt_diferencia_pre ht_dhdt_medido_pre parametros_acum_pre
    clear file_temp_puntual Tc_dem
end
%matriz_de_comprobacion=[codigos_simulados dif_vols_mod_totales dif_area_obs_totales a_obs_ac'];
% Datos para análisis
save(strcat(output_link,'/stat_test.mat'),'std_data_for_map');
save(strcat(output_link,'/No_modelado.mat'),'glaciares_no_modelados');
%% 10.7 Cierre de script
tEnd = toc(tStart)
