clear all
tStart = tic;
grayColor = [.7 .7 .7];

% Link datos georeferenciales para cada glaciar
link_gl = '/directorio/info';

% Link a resumen de modelación.
link_db='/directorio/modelo';

% directorio archivos de forzantes meteorológicas
link_meteo=strcat('/directorio/forzantes');

%Guardar archivos

model='/nombre_modelo'; % Nombre de Modelo

main='/directorio_de_salida'; % Directorio Principal de Salida de archivos

mkdir(strcat(main,model)) % Creación de directorio

output_link=strcat(main,model);

mkdir(strcat(output_link,'/Qg_y_params')); % Creación de directorio

link_save=strcat(output_link,'/Qg_y_params'); % Directorio de salida de modelo para calibración

mkdir(strcat(output_link,'/csv')); % Creación de directorio

link_publico=strcat(output_link,'/csv'); % Directorio de archivos de salida

% Se cargan vectores de latitud y longitud para referenciar
% coordenadas <REASIGNAR DIRECTORIO>.
load /directorio/geodesico_referencial/lat

load /directorio/geodesico_referencial/lon


load(strcat(link_db,'/stat_test','.mat'));

std_data_for_map(std_data_for_map(:,1)==0,:)=[];

codigos=std_data_for_map(:,1);  % Carga de codigo glaciar COD_GLA

gl_modelados=length(codigos); 

%Parametrización Física
factor_area = 0.5:0.1:1.5; % Factor de calibración de área

b=1.3598;            % Factor de ajuste Volume Area Scaling

c=0.027;             % Factor de ajuste Volume Area Scaling

% Constante de Stefann-Boltzmann
sigma = 5.67*10^(-8); % [ W / (m2 * K4)]

% Emisividad de la nieve. Rango: 0.97 a 1.00
snow_emissivity = 0.99;

% Emisividad del aire.
air_emissivity = 1;

% Temperatura de la superficie del hielo y manto
T_snow = -1;
% albedo máximo y mínimo que de la nieve y el glaciar
alb_max_snow = 0.90;

alb_min_snow = 0.60;

alb_ice = 0.40;       % Albedo del hielo

% Temperatura de separación lluvia/nieve
To= 0; % [C] Temperatura para umbral de derretimiento

r = 0.10; % [1/d]

glaciares_no_modelados = []; %Creación de matriz vacia para modelados fallidos

% Parametrización para balance de energía
k_o=0.4;                    %constante de von karmann

z= 2.0;                     %Altura de forzante

z_o= 0.0005;                % Parametro de rugosidad (extraido de Cuffey and paterson 2010)

C=((k_o^2)/log(z/z_o)^2);   % C*  Coeficiente de transferencia efectiva

t_s=0;                      %Temperatura superficial °C

Avps=611;                   %¨Presion de vapor en la superificie [Pascal]

%Vectores de tiempo
tiempo_day = datenum(1981,01,01) : 1 : datenum(2020,04,30); %Periodo acotado para modelacion

tiempo_day = tiempo_day';

yr_mod=1981; %Año de inicio de simulación

tiempo_plots = datetime(tiempo_day,'ConvertFrom','datenum');

tiempo_tablas = datevec(tiempo_day);

yr=tiempo_tablas(:,1);

mth=tiempo_tablas(:,2);

day=tiempo_tablas(:,3);

l_inicio=1;

l_final=gl_modelados;

% Inicio de modelado

for i=l_inicio:l_final

    %Carga de datos
    load(strcat(link_gl,'/GLA_CL',num2str(codigos(i)),'.mat'))

    fila = find(abs(lat-info(2)) == min(abs(lat-info(2))));

    col = find(abs(lon-info(3)) == min(abs(lon-info(3))));

    lat1=lat(fila);    %Latitud

    lon1=lon(col);     %Longitud

    area_inv=std_data_for_map(i,6); %Área inventario

    area_yr=std_data_for_map(i,8);  %Año de inventario

    area_yr_take=std_data_for_map(i,5); %Año de inventario ajustado

    area_cal=area_inv*factor_area;

    try

    Tc=info(23);

    catch

    Tc=0;

    end
    % Indices para validacion
    index_ainv=find(tiempo_day==datenum(area_yr,01,01));

    % Carga de forzantes
    if exist(strcat(link_meteo,'/data_',num2str(lat(fila)),'_',num2str(lon(col))),'file') > 0 || exist(strcat(link_meteo,'/data_CL',num2str(info(1)),'.txt'),'file')>0 %Dato

        % Si el glaciar tiene un área mayor a 25 km^2 se carga temperatura de promedio de varios puntos
        if info(:,7)/10^6>=25
            file_temp_poligonal=readmatrix(strcat(link_meteo,'/data_CL',num2str(info(1)),'.txt'));

            temperatura=file_temp_poligonal(:,7);

            pr=file_temp_poligonal(:,2);

            ws=file_temp_poligonal(:,3);

            pa=file_temp_poligonal(:,4);

            sw=file_temp_poligonal(:,5);

            lw=file_temp_poligonal(:,6);

          % Corrección mediante gradiente de temperatura -6.49 K/km ISA (international standar atmosphere) 0-11000 msnm

            t2m=file_temp_poligonal(:,7)+Tc;

            rh=file_temp_poligonal(:,8);

            % Revisión de input de temperatura
            if isnan(temperatura(1))==false
                % Si el archivo no contiene información se carga el dato de temperatura a partide un solo punto.

            elseif exist('temperatura')==false
                file_temp_puntual=dlmread(strcat(link_meteo,'/data_',num2str(lat(fila)),'_',num2str(lon(col))),'\t');

                pr=file_temp_puntual(:,1);

                ws=file_temp_puntual(:,2);

                pa=file_temp_puntual(:,3);

                sw=file_temp_puntual(:,4);

                lw=file_temp_puntual(:,5);


                t2m=file_temp_puntual(:,6)+Tc;

                rh=file_temp_puntual(:,7);



            end % En de if (231)
        else
            file_temp_puntual=load(strcat(link_meteo,'/data_',num2str(lat(fila)),'_',num2str(lon(col))),'\t');

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

    end
    % Estructura de archivo
    % (1) info(1) (2) 1 (3)lon1 (4)lat1 (5)info(6) (6)a_obs (7)v_obs (8)yr_mod ...
    % (9)parametros_gl(1,7) (10)parametros_gl(1,8) (11)bn_dhdt_medido
    % (12)bn_dhdt_modelado (13)bn_dhdt_diferencia (14)bn_dhdt_diferencia_r
    % (15)ht_dhdt_medido (16)ht_dhdt_modelado (17)ht_dhdt_diferencia
    % (18)ht_dhdt_diferencia_r (19)Tc (20)Ts];

    %% Comienza modelamiento
    if std_data_for_map(i,2)==1  %if zonal

        zona="NORTE";

        c1=std_data_for_map(i,9);

        c0=std_data_for_map(i,10);

        ts=std_data_for_map(i,20);

        sw_f=std_data_for_map(i,21);

        tc=0;

        pairs_zn=[c1,c0,sw_f];

        vapourpressuresaturation=(6.116441*10.^((7.591386.*t2m)./(t2m+240.7263))).*100;

        vapourpressure=(rh.*vapourpressuresaturation);

        % Input para modelación de derretimiento

        % [Version completa] x_diario=[tiempo_day pr t2m rh pa ws vapourpressure sw lw];
        x_diario=[tiempo_day pr t2m rh pa ws vapourpressure sw lw]; % Versión acotada

        for qq=1:length(area_cal)

            if info(1,13) == 1 % verificar si el glaciar es rocoso

                [parametros, flujos] = balance_ZN(x_diario, To, ts, area_cal(qq), c, b, pairs_zn, alb_ice, 0.2, C, Avps,t_s,yr_mod ,info(22));

            else

                [parametros, flujos] = balance_ZN(x_diario, To, ts, area_cal(qq), c, b, pairs_zn, alb_ice, 1, C, Avps,t_s,yr_mod ,info(22));

            end

            flujos_acum_pre(:,:,qq)=flujos;         %Resultados previos para cada Ts posible

            parametros_acum_pre(:,:,qq)=parametros; %Resultados previos para cada Ts posible

        end

    elseif std_data_for_map(i,2)==2  %if zonal

        zona="SUR";

        DDI=std_data_for_map(i,9);

        DDS=std_data_for_map(i,10);

        ts=std_data_for_map(i,20);

        tc=0;

        pairs=[DDI,DDS];



        x_diario=[tiempo_day pr t2m];

        % Comienza modelado
        for qq=1:length(area_cal)

            if info(1,13) == 1 % verificar si el glaciar es rocoso

                [parametros, flujos] = balance_ZA(x_diario, To, ts,  area_cal(qq), c, b, pairs, alb_ice, 0.2, yr_mod, info(22));

            else

                [parametros, flujos] = balance_ZA(x_diario, To, ts, area_cal(qq), c, b, pairs, alb_ice, 1, yr_mod, info(22));

            end

            flujos_acum_pre(:,:,qq)=flujos;         %Resultados previos para cada Ts posible

            parametros_acum_pre(:,:,qq)=parametros; %Resultados previos para cada Ts posible

        end


    end %if zonal
    
    % Busqueda de área modelada más proxima al área inventariada para
    % el glaciar modelado.
    for s=1:size(flujos_acum_pre,3)

    dif_area(:,s)=area_inv-flujos_acum_pre(index_ainv,2,s);

    dif_abs(:,s)=abs(area_inv-flujos_acum_pre(index_ainv,2,s));

    end

    best_dif=min(dif_abs);

    best_dif_index=find(dif_abs==best_dif);

    if length(best_dif_index)>1

        best_dif_index=best_dif_index(1);

    else

        %nada

    end

   % Guardando mejor salida del modelo.

    flujos_gl=flujos_acum_pre(:,:,best_dif_index);

    parametros_gl=parametros_acum_pre(:,:,best_dif_index);

    parametros_gl=[parametros_gl,best_dif];

%% Cracion de tablas
   
   % Balance de masa
   dhdt=flujos_gl(:,3)./flujos_gl(:,2);

   %Guardando resultado para DGA
   if std_data_for_map(i,2)==1
                    Cabeceras_diarias = {'Año';'Mes';'Día';'SW [W/m^2]';'Temperatura [°C]';'Precipitación [mm]';'Precipitación sólida[mm]';'Área[km2]';'Volumen[km3]';'SWE [mm]';'Derretimiento nieve[mm]';'Derretimiento hielo [mm]';'Caudal Nival[m3/s]';'Caudal Hielo[m3/s]';'Caudal Total[m3/s]';'Sublimación [mm]';'dh/dt [m/yr]'};

                    tabla_diaria=table(yr,mth,day,t2m,pr,flujos_gl(:,13),flujos_gl(:,2),flujos_gl(:,3),flujos_gl(:,8),flujos_gl(:,9),flujos_gl(:,7),flujos_gl(:,14),flujos_gl(:,12),flujos_gl(:,15),flujos_gl(:,end),dhdt,'VariableNames',Cabeceras_diarias);

                   % Cabeceras_info={'Codigo Glaciar';'Latitud [°]';'Longitud [°]';'Altitud [msnm]';'Área inventario [km2]';'Zona glaciológica';'Rocoso [boolean]';'dh/dt[m/año]';'Suma de errores porc. [%]';'C1';'C0'};
                   % tabla_informacion_glaciar=table(info(1),info(2),info(3),info(1,12),(info(:,7)/10^6),string(zona),info(1,13),info(1,14),puntaje(puntaje_menor),parametros_gl(1,7),parametros_gl(1,8),'VariableNames',Cabeceras_info);

   elseif std_data_for_map(i,2)==2
                    Cabeceras_diarias = {'Año';'Mes';'Día';'SW [W/m^2]';'Temperatura [°C]';'Precipitación [mm]';'Precipitación sólida[mm]';'Área[km2]';'Volumen[km3]';'SWE [mm]';'Derretimiento nieve[mm]';'Derretimiento hielo [mm]';'Caudal Nival[m3/s]';'Caudal Hielo[m3/s]';'Caudal Total[m3/s]';'Sublimación [mm]';'dh/dt [m/yr]'};

                    tabla_diaria=table(yr,mth,day,t2m,pr,flujos_gl(:,13),flujos_gl(:,2),flujos_gl(:,3),flujos_gl(:,8),flujos_gl(:,9),flujos_gl(:,7),flujos_gl(:,14),flujos_gl(:,12),flujos_gl(:,15),NaN,dhdt,'VariableNames',Cabeceras_diarias);

                   % Cabeceras_info={'Codigo Glaciar';'Latitud [°]';'Longitud [°]';'Altitud [msnm]';'Área inventario [km2]';'Zona glaciológica';'Rocoso [boolean]';'dh/dt[m/año]';'Suma de errores porc. [%]';'DDI';'DDS'};
   end

   % Escrituta de archivos xlsx para DGA
   writetable(tabla_diaria,strcat(link_publico,'/CL',num2str(codigos(i)),'.csv'),'Delimiter',',');
                    % writetable(tabla_informacion_glaciar,strcat(link_rpublico,'/CL_',string(Glaciares_prueba.nombre(i)),'_',string(zona),'_',prefijo,'.xlsx'),'Sheet',2);


   % Resumen

   if lat1>-36 %1 Norte

       resumen(i,:)=[info(1) 1 lon1 lat1 info(6) area_inv parametros_gl(1,7) parametros_gl(1,8) parametros_gl(:,end) parametros_gl(1,1)  best_dif];

   else lat1<=-36; % 2 Sur
       resumen(i,:)=[info(1) 2 lon1 lat1 info(6) area_inv parametros_gl(1,7) parametros_gl(1,8) NaN parametros_gl(1,1) best_dif];

   end

    % Guardando salida para balance hidrico
    save(strcat(link_save,'/Q_y_params_GL',num2str(codigos(i))), 'parametros_gl' , 'flujos_gl')

    clear area_yr_take best_dif best_dif_index c0 c1 dif_abs dif_area fila ...
        file_temp_puntual flujos flujos_acum_pre flujos_gl index_ainv ...
        info lat1 lon1 lw pairs_zn pairs parametros_acum_pre parametros_gl ...
        pr qq s sw t2m tc ts vapourpressure vapourpressuresaturation ws ...
        x_diario z zona tabla_diaria Cabeceras_diarias sw_f
    disp(i+" de "+gl_modelados)
end

tEnd = toc(tStart)

save(strcat(output_link,'/resumen.mat'),'resumen');
