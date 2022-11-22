clear all
tStart = tic;
grayColor = [.7 .7 .7];
% Link datos georeferenciales para cada glaciar

link_gl = '/home/duiliof/Desktop/BH/input/info_HBH_filtrado';
% Link a resumen de modelación.
link_db='/media/duiliof/4T_HDD1/RESPALDO_BH5/OUTPUTS/BASE_OUTPUT/V25_TcDEM_hnmean';
% Link meteo
link_meteo=strcat('/home/duiliof/Desktop/BH/input/GCM/IPSL-CM5A-LR');
%Guardar archivos
model='/IPSL-CM5A-LR';
main='/home/duiliof/Desktop/BH/resultados';
mkdir(strcat(main,model))
output_link=strcat(main,model);
mkdir(strcat(output_link,'/Qg_y_params'));
link_save=strcat(output_link,'/Qg_y_params');
mkdir(strcat(output_link,'/csv'));
link_publico=strcat(output_link,'/csv');
S=shaperead("/home/duiliof/Desktop/BH/input/Inventario/IPG2014/inv_ful_2014_EPSG4326_BRAUN_HUGONNET_filtrado_vecinos");
%S=shaperead("/home/duiliof/Desktop/BH/input/Inventario/IPG2014/inventario_2014_bandas_hugonnet_vecinos.shp");

% Se cargan vectores de latitud y longitud para referenciar
% coordenadas <REASIGNAR DIRECTORIO>.
load /home/duiliof/Desktop/BH/input/geodesico_referencial/lat
load /home/duiliof/Desktop/BH/input/geodesico_referencial/lon


load(strcat(link_db,'/stat_test','.mat'));
std_data_for_map(std_data_for_map(:,1)==0,:)=[];
codigos=std_data_for_map(:,1);
gl_modelados=length(S);

%Parametrización Física
factor_area = 0.5:0.1:1.5;
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
Tc= 0;
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
tiempo_day = datenum(1981,01,01) : 1 : datenum(2100,12,31); %Periodo acotado para modelacion
tiempo_day = tiempo_day';
yr_mod=1981; %Año de inicio de simulación
tiempo_plots = datetime(tiempo_day,'ConvertFrom','datenum');
tiempo_tablas = datevec(tiempo_day);
yr=tiempo_tablas(:,1);
mth=tiempo_tablas(:,2);
day=tiempo_tablas(:,3);
l_inicio=1;
l_final=length(S);


% Consulta por glaciar
%prompt=["GLACIARES PARA PRUEBA";"Todos           (0)";"Parinacota      (1)";"Tronquitos      (2)";"San Francisco   (3)";"Bello           (4)";"Mocho           (5)";"Tyndall         (6)";"Exploradores    (7)"];
%disp(prompt)
%glaciar_prueba=input('Por favor seleccione el/los glaciares para modelar:');
%zona_var=input('Habilitar zonas glaciológicas automática 0(NO) o 1(SI):');
% if zona_var==0
%     zona_forzada=input('Por favor ingrese zona que desea forzar entre comillas(ZN/ZC/ZS/ZA):');
% else
%     zona_forzada=NaN;
% end



for i=l_inicio:l_final
%     for i=35
    %Carga de datos
    %load(strcat(link_gl,'/GLA_CL',num2str(Glaciares_prueba.codigo(i)),'.mat'))
    load(strcat(link_gl,'/GLA_',S(i).COD_GLA,'.mat'))
    index_db=find(S(i).vecinos==string(strcat('CL',num2str(std_data_for_map(:,1)))));
    fila_2000 = find(abs(lat-info(2)) == min(abs(lat-info(2))));
    col = find(abs(lon-info(3)) == min(abs(lon-info(3))));
    lat1=lat(fila_2000);    %Latitud
    lon1=lon(col);          %Longitud
    area_inv=(info(:,7)/10^6); %Área inventario
    area_yr=info(6);  %Año de inventario
    area_cal=area_inv*factor_area;
    %Tc load
    try 
    Tc=info(23);
    catch
    Tc=0; 
    end
    %Indices para validacion
    index_ainv=find(tiempo_day==datenum(area_yr,01,01));
    
    %Carga de forzantes
    try
    if exist(strcat(link_meteo,'/data_',num2str(lat(fila_2000)),'_',num2str(lon(col))),'file') > 0 || exist(strcat(link_meteo,'/data_CL',num2str(info(1)),'.txt'),'file')>0 %Dato
        
        % Si el glaciar tiene un área mayor a 25 km^2 se carga temperatura de promedio de varios puntos
        if info(:,7)/10^6>=25
            file_temp_poligonal=readmatrix(strcat(link_meteo,'/data_CL',num2str(info(1)),'.txt'));
            temperatura=file_temp_poligonal(:,4);
            pr=file_temp_poligonal(:,2);
            sw=file_temp_poligonal(:,3);
            t2m=file_temp_poligonal(:,4)+Tc;
            
            
            % Revisión de input de temperatura
            if isnan(temperatura(1))==false
                % 7.4. Si el archivo no contiene información se carga el dato de temperatura a partide un solo punto.
            elseif exist('temperatura')==false
                file_temp_puntual=importdata(strcat(link_meteo,'/data_',num2str(lat(fila_2000)),'_',num2str(lon(col))),'\t',1);
                pr=file_temp_puntual.data(:,1);
                sw=file_temp_puntual.data(:,3);
                t2m=file_temp_puntual.data(:,2)+Tc;
                % 7.5. Corrección mediante gradiente de temperatura -6.49 K/km ISA (international standar atmosphere) 0-11000 msnm
                
            end % En de if (231)
        else
            file_temp_puntual=importdata(strcat(link_meteo,'/data_',num2str(lat(fila_2000)),'_',num2str(lon(col))),'\t',1);
            pr=file_temp_puntual.data(:,1);
            sw=file_temp_puntual.data(:,3);
            t2m=file_temp_puntual.data(:,2)+Tc;
            
        end
        
    end
    %Estructura de archivo
    % (1) info(1) (2) 1 (3)lon1 (4)lat1 (5)info(6) (6)a_obs (7)v_obs (8)yr_mod ...
    % (9)parametros_gl(1,7) (10)parametros_gl(1,8) (11)bn_dhdt_medido
    % (12)bn_dhdt_modelado (13)bn_dhdt_diferencia (14)bn_dhdt_diferencia_r
    % (15)ht_dhdt_medido (16)ht_dhdt_modelado (17)ht_dhdt_diferencia
    % (18)ht_dhdt_diferencia_r (19)Tc (20)Ts];
    
    %% Comienza modelamiento
    if std_data_for_map(index_db,2)==1  %if zonal
        zona="NORTE";
        c1=std_data_for_map(index_db,9);
        c0=std_data_for_map(index_db,10);
        c2=1;
        ts=std_data_for_map(index_db,20);
        tc=0;
        pairs_zn=[c1,c0,c2];
               
        % Input para modelación de derretimiento
        % [Version completa] x_diario=[tiempo_day pr t2m rh pa ws vapourpressure sw lw];
        x_diario=[tiempo_day pr t2m sw]; % Versión acotada
        
        for qq=1:length(area_cal)
            if info(1,13) == 1 % verificar si el glaciar es rocoso
                [parametros, flujos] = balance_ZN_tiempo_acotado(x_diario, To, ts, area_cal(qq), c, b, pairs_zn, alb_ice, 0.2, C, Avps,t_s,yr_mod ,info(22));
            else
                [parametros, flujos] = balance_ZN_tiempo_acotado(x_diario, To, ts, area_cal(qq), c, b, pairs_zn, alb_ice, 1, C, Avps,t_s,yr_mod ,info(22));
            end
            flujos_acum_pre(:,:,qq)=flujos;         %Resultados previos para cada Ts posible
            parametros_acum_pre(:,:,qq)=parametros; %Resultados previos para cada Ts posible
        end
        
    elseif std_data_for_map(index_db,2)==2  %if zonal
        zona="SUR";
        DDI=std_data_for_map(index_db,9);
        DDS=std_data_for_map(index_db,10);
        ts=std_data_for_map(index_db,20);
        tc=0;
        pairs=[DDI,DDS];
        
        %Modelamiento
        x_diario=[tiempo_day pr t2m];
        
        for qq=1:length(area_cal)
            if info(1,13) == 1 % verificar si el glaciar es rocoso
                [parametros, flujos] = balance_ZA_tiempo_acotado(x_diario, To, ts,  area_cal(qq), c, b, pairs, alb_ice, 0.2, yr_mod, info(22));
            else
                [parametros, flujos] = balance_ZA_tiempo_acotado(x_diario, To, ts, area_cal(qq), c, b, pairs, alb_ice, 1, yr_mod, info(22));
            end
            flujos_acum_pre(:,:,qq)=flujos;         %Resultados previos para cada Ts posible
            parametros_acum_pre(:,:,qq)=parametros; %Resultados previos para cada Ts posible
        end
        
    end %if zonal
    
   
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
    
%     % Guardando salidas del modelo.
    flujos_gl=flujos_acum_pre(:,:,best_dif_index);
    parametros_gl=parametros_acum_pre(:,:,best_dif_index);
    parametros_gl=[parametros_gl,best_dif];
%     figure
%     for z=1:size(flujos_acum_pre,3)
%         plot(tiempo_plots,flujos_acum_pre(:,2,z),'Color', grayColor)
%         hold on
%     end
%        plot(tiempo_plots,flujos_acum_pre(:,2,best_dif_index),'Color','red');
%        scatter(tiempo_plots(index_ainv),area_inv,100,'o','red','filled');
%        grid on
%        xlabel('Año');
%        ylabel('Área [km^2]');
%        hold off

%% Cracion de tablas 
   
   %Guardando resultado para DGA
   %% Cracion de tablas 
   dhdt=flujos_gl(:,3)./flujos_gl(:,2);
   %Guardando resultado para DGA
   if std_data_for_map(i,2)==1
                    Cabeceras_diarias = {'Año';'Mes';'Día';'SW [W/m^2]';'Temperatura [°C]';'Precipitación [mm]';'Precipitación sólida[mm]';'Área[km2]';'Volumen[km3]';'SWE [mm]';'Derretimiento nieve[mm]';'Derretimiento hielo [mm]';'Caudal Nival[m3/s]';'Caudal Hielo[m3/s]';'Caudal Total[m3/s]';'Sublimación [mm]';'dh/dt [m/yr]'};
                    tabla_diaria=table(yr,mth,day,sw,t2m,pr,flujos_gl(:,13),flujos_gl(:,2),flujos_gl(:,3),flujos_gl(:,8),flujos_gl(:,9),flujos_gl(:,7),flujos_gl(:,14),flujos_gl(:,12),flujos_gl(:,15),flujos_gl(:,end),dhdt,'VariableNames',Cabeceras_diarias);
                   % Cabeceras_info={'Codigo Glaciar';'Latitud [°]';'Longitud [°]';'Altitud [msnm]';'Área inventario [km2]';'Zona glaciológica';'Rocoso [boolean]';'dh/dt[m/año]';'Suma de errores porc. [%]';'C1';'C0'};
                   % tabla_informacion_glaciar=table(info(1),info(2),info(3),info(1,12),(info(:,7)/10^6),string(zona),info(1,13),info(1,14),puntaje(puntaje_menor),parametros_gl(1,7),parametros_gl(1,8),'VariableNames',Cabeceras_info);
   elseif std_data_for_map(i,2)==2
                    Cabeceras_diarias = {'Año';'Mes';'Día';'SW [W/m^2]';'Temperatura [°C]';'Precipitación [mm]';'Precipitación sólida[mm]';'Área[km2]';'Volumen[km3]';'SWE [mm]';'Derretimiento nieve[mm]';'Derretimiento hielo [mm]';'Caudal Nival[m3/s]';'Caudal Hielo[m3/s]';'Caudal Total[m3/s]';'Sublimación [mm]';'dh/dt [m/yr]'};
                    tabla_diaria=table(yr,mth,day,sw,t2m,pr,flujos_gl(:,13),flujos_gl(:,2),flujos_gl(:,3),flujos_gl(:,8),flujos_gl(:,9),flujos_gl(:,7),flujos_gl(:,14),flujos_gl(:,12),flujos_gl(:,15),nan(length(yr),1),dhdt,'VariableNames',Cabeceras_diarias);
                   % Cabeceras_info={'Codigo Glaciar';'Latitud [°]';'Longitud [°]';'Altitud [msnm]';'Área inventario [km2]';'Zona glaciológica';'Rocoso [boolean]';'dh/dt[m/año]';'Suma de errores porc. [%]';'DDI';'DDS'};
   end
   
   % Escrituta de archivos xlsx para DGA
   writetable(tabla_diaria,strcat(link_publico,'/',S(i).COD_GLA,'.csv'),'Delimiter',',');
                    % writetable(tabla_informacion_glaciar,strcat(link_rpublico,'/CL_',string(Glaciares_prueba.nombre(i)),'_',string(zona),'_',prefijo,'.xlsx'),'Sheet',2);
   %% Balance de Masa
periodo_de_validacion=(datenum(1981,01,01,00,00,00):datenum(0000,00,01,00,00,00):datenum(2020,04,30,00,00,00))';
MB=table();
MB.area=flujos_gl(:,2);
MB.volumen=flujos_gl(:,3);
yr_frac=(periodo_de_validacion(end)-periodo_de_validacion(1))/365.25;
dhdt_modelado=(MB.volumen(end)-MB.volumen(1))*1000/MB.area(1)/yr_frac; 
                    
   % Resumen 
   
  if lat1>-36 %1 Norte
       resumen(i,:)=[info(1) 1 lon1 lat1 info(6) area_inv parametros_gl(1,7) parametros_gl(1,8) parametros_gl(:,end) parametros_gl(1,1)  best_dif Tc dhdt_modelado];
   else lat1<=-36; % 2 Sur
       resumen(i,:)=[info(1) 2 lon1 lat1 info(6) area_inv parametros_gl(1,7) parametros_gl(1,8) NaN parametros_gl(1,1) best_dif Tc dhdt_modelado];
   end
    % Guardando salida para balance hidrico
    save(strcat(link_save,'/Q_y_params_GL',S(i).COD_GLA), 'parametros_gl' , 'flujos_gl')
    clear area_yr_take best_dif best_dif_index c0 c1 dif_abs dif_area fila_2000 ...
        file_temp_puntual flujos flujos_acum_pre flujos_gl index_ainv ...
        info lat1 lon1 lw pairs_zn pairs parametros_acum_pre parametros_gl ...
        pr qq s sw t2m tc ts vapourpressure vapourpressuresaturation ws ...
        x_diario z zona tabla_diaria Cabeceras_diarias 
    clc
    disp(i+" de "+gl_modelados)
    
    catch
        no_data(i,1)=string([S(i).COD_GLA]);
    end 
end
tEnd = toc(tStart)
save(strcat(output_link,'/resumen.mat'),'resumen');
