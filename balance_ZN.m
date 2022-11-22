function [params, flujos] = balance_ZN(forz, To, Ts, Ai, cc, bb, pairs_zn, a_ice, rocoso, C, Avps, t_s, yr_mod, calving)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parámetros
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forz: contiene la fecha, forziación de onda larga neta, onda corta
% incidente, precipitación diaria y temperatura [J/m2]
% To: Temperatura que separa la lluvia de la nieve [ºC]
% r : tasa de decaimiento del Albedo [1/dias]
% Ai: Área inicial del glaciar [km2]
% cc: Intercambio volumen glaciar
% bb: Intercambio volumen glaciar
% a_ice: albedo del glaciar
% Rf: factor que modifica la forzante radiativa
tiempo_day2 = forz(1,1) : 1 : forz(end, 1);
tiempo_day_output=tiempo_day2;
tiempo_day2=datevec(tiempo_day2);
sw_f=pairs_zn(1,3);
% fecha en la cual se realizará el cambio de Área y Volumen
actualizar_V = [3, 15]; % [mes dia]
%% Forzantes
pr=forz(:,2);   %Precipitación [mm]
t2m=forz(:,3);  %temepratura [C]
t2m_K=forz(:,3)+273.15;  %temepratura [K]
sw=forz(:,4);   %forziación de onda corta [W/m2]

%%  Parámetros físicos
r = 0.12; % [1/d]
%r = 0.10; % [1/d]
rho_w = 1000; % [kg/m3]
rho_ice=917; %[km/m3]
rho_snow=400; %[km/m3]
Lf = 334000; %[J/kg] Calor latente de fusi�n
Ls = 2830000; %[J/kg] Calor latente de sublimación
B = 0.95; % Rango 0.8 a 1.1
Fice = 1; % Proporci�n del �rea superficial con hielo.
albedo_ice = a_ice; % albedo que tiene el hielo una vez que no hay nieve
k_ice=3.5; % factor que retrasa el derretimiento del hielo
k_snow=12.5; % factor que retrasa el derretimeinto de la nieve
%albedo_snow_min = 0.4;
%albedo_snow_max = 0.85;
albedo_snow_min = 0.6;
albedo_snow_max = 0.9;
reduccion_forz_rain = 0.3; % reducci�n de la forziaci�n en d�as con lluvia DeWalle&Rango
lw_out=-315.6; % Lw out provisional asumiendo tempereratura constante de 273.15K

% Caso de calving
if isnan(calving)==true
    calving=0;
else 
    calving=calving;
end 

%% Crear vectores con variables de estado
subl= zeros(length(forz(:,1)), 1);
area_gl = zeros(length(forz(:,1)), 1);
volumen = zeros(length(forz(:,1)), 1);
Sfree = zeros(length(forz(:,1)), 1);
albedo_snow = zeros(length(forz(:,1)), 1);
Mpot_snow = zeros(length(forz(:,1)), 1);
Mpot_ice = zeros(length(forz(:,1)), 1);
SUBLpot_snow = zeros(length(forz(:,1)), 1);
SUBLpot_ice = zeros(length(forz(:,1)), 1);
SWE = zeros(length(forz(:,1)), 1);
M_snow = zeros(length(forz(:,1)), 1);
EE= zeros(length(forz(:,1)), 1);
Qsnow = zeros(length(forz(:,1)), 1);
Qice = zeros(length(forz(:,1)), 1);
SWnet_snow=zeros(length(forz(:,1)), 1);
SWnet_ice=zeros(length(forz(:,1)), 1);
Pliq = zeros(length(forz(:,1)), 1);
Psnow = zeros(length(forz(:,1)), 1);
Em_snow = zeros(length(forz(:,1)), 1);
Em_ice = zeros(length(forz(:,1)), 1);
Es_snow = zeros(length(forz(:,1)), 1);
Es_ice = zeros(length(forz(:,1)), 1);
Sublimacion=zeros(length(forz(:,1)), 1);
% Iniciar variables de estado
area_gl(1) = Ai; % Se asigna el �rea inicial en km2
volumen(1) = cc * (Ai)^bb; % Se obtiene el volume inicial en km3
Sfree(1) = 0.5; % Se asigna el valor de Sfree inicial
albedo_snow(1) = 0.5; % albedo inicial
C1=pairs_zn(1,1); %Factor de calibración
C0=pairs_zn(1,2); %Factor de calibración

% parámetros
% días desde que cay� la ultima nevada
n = 0;


%% Análisis diario

Vice = 0;

ver = [];

conteo_agnos_SWE = 0;

%SWE calculo inicial

for xx = 2:1825
    
    if t2m(xx) > Ts
        Pliq(xx) = pr(xx); % se asigna lo que cae como lluvia
    else
        Psnow(xx) = pr(xx); % se asigna lo que cae como nieve
        
        if pr(xx) > 5 % 5 [mm]/d
            n = 0; % se resetea el conteo de nieve
        else
            n = n+1;
        end
    end
    
    % Decaimiento del albedo

 %% 2. se calcula el albedo de la nieve
  
        albedo_snow(xx) = albedo_snow_min + (albedo_snow_max - albedo_snow_min) * exp(-r * n);    
%     
    if pr(xx) > 5 % 5 mm/d
        % se corrige como día con precipitación
        SWnet_snow(xx) = (1-albedo_snow(xx)) * sw(xx) * reduccion_forz_rain ; % [J/m2]
        SWnet_ice(xx) = (1-albedo_ice) * sw(xx) * reduccion_forz_rain ; % [J/m2]
    else
        SWnet_snow(xx) = (1-albedo_snow(xx)) * sw(xx) ; % [J/m2]
        SWnet_ice(xx) = (1-albedo_ice) * sw(xx); % [J/m2]
    end


   Em_ice(xx)=(SWnet_ice(xx)*sw_f+(t2m(xx)*C1)+C0)*(1000*3600*24/(Lf*rho_w));
   Em_snow(xx)=(SWnet_snow(xx)*sw_f+(t2m(xx)*C1)+C0)*(1000*3600*24/(Lf*rho_w));
   
       
    %% 4. Se calculan los derretimientos potenciales
    
    if t2m(xx) > To
        Mpot_snow(xx) = Em_snow(xx) ; % [mm]
        Mpot_ice(xx) = rocoso * Sfree(xx-1) * Em_ice(xx) ;% [mm]
        
    else
        Mpot_snow(xx) = 0 ;
        Mpot_ice(xx) = 0 ;
    end
    %
    if Mpot_snow(xx) < 0
        Mpot_snow(xx) = 0;
    end
    if Mpot_ice(xx) < 0
        Mpot_ice(xx) = 0;
    end
    
    
    %% 5. Se obtiene el monto de SWE y el derretimiento
    
    SWE(xx) = SWE(xx-1) + Psnow(xx);
    M_snow(xx) = min(SWE(xx), Mpot_snow(xx));
    
    %% 6. Obtener los caudales de la nieve y el hielo
    
    Qsnow(xx) = Qsnow(xx-1) * exp(-1/k_snow) + (Pliq(xx) + M_snow(xx)) * (1-exp(-1/k_snow)); % [mm]
    
    Qice(xx) = Qice(xx-1) * exp(-1/k_ice) + (Pliq(xx) * Sfree(xx) + Mpot_ice(xx)) * (1 - exp(-1/k_ice)); % [mm]
    
    %% 7. Actualizar el monto de SWE sacando lo que se derriti� y Sfree
    
    SWE(xx) = SWE(xx) - M_snow(xx);
    
    if SWE(xx) < 0
        SWE(xx) = 0;
    end
    % Actualizar Sfree
    if Psnow(xx) > 0
        Sfree(xx) = 0;
    elseif SWE(xx) == 0
        Sfree(xx) = 1;
    elseif 1 - (SWE(xx) - M_snow(xx)) / SWE(xx) >= 0 & 1 - (SWE(xx) - M_snow(xx)) / SWE(xx) <= 1
        Sfree(xx) = 1 - (SWE(xx) - M_snow(xx)) / SWE(xx) ;
    else
        Sfree(xx) = 1;
    end
    
    %% 8. Actualización del volumen

    
    if tiempo_day2(xx,2) == actualizar_V(1) & tiempo_day2(xx,3) == actualizar_V(2)
        
        % Se incorpora el SWE si ha duforzo muchos años
        if SWE(xx) > 0
            conteo_agnos_SWE = conteo_agnos_SWE + 1;
        else
            conteo_agnos_SWE = 0;
        end
        
        if conteo_agnos_SWE == 1
            % Se vuelve a resetear el contador
            conteo_agnos_SWE = 0;
            
            volumen(xx) = volumen(xx-1) - Vice + calving + (1-Sfree(xx))*SWE(xx)*area_gl(xx-1)/10^6;
            % todo el SWE pasar a ser hielo
            SWE(xx) = 0;
        else
            volumen(xx) = volumen(xx-1) - Vice + calving;
            
        end
        
        ver = [ver Vice];
        
        
        area_gl(xx) = (volumen(xx) / cc ) ^ (1/bb);
        if volumen(xx) < 0
            volumen(xx) = 0;
            area_gl(xx) = 0;
        end
        Vice = 0;
    else
        area_gl(xx) = area_gl(xx-1);
        volumen(xx) = volumen(xx-1);
        Vice = Vice + ( area_gl(xx)*1000*1000 * Mpot_ice(xx)/1000 ) / 10^9;
    end
    
    
end
for zz=1:10
% [Corrida completa]    init_date=find(tiempo_day_output==datenum(1981+zz-1,12,01));
% [Corrida completa]    final_date=find(tiempo_day_output==datenum(1981+zz,03,14));
% [Corrida completa]    mean_swe_year(zz,1)=mean(SWE(init_date:final_date));
    init_date=find(tiempo_day_output==datenum(yr_mod+zz-1,03,16));
    final_date=find(tiempo_day_output==datenum(yr_mod+zz,12,31));
    mean_swe_year(zz,1)=mean(SWE(init_date:final_date));
end

clear volumen area_gl Vice ver SWE conteo_agnos_SWE Sfree Qice Qsnow Mpot_snow Mpot_ice Em_ice Em_snow Pliq Psnow SWnet_ice SWnet_snow

% Empieza el modelo
% Reload variables
subl= zeros(length(forz(:,1)), 1);
area_gl = zeros(length(forz(:,1)), 1);
volumen = zeros(length(forz(:,1)), 1);
Sfree = zeros(length(forz(:,1)), 1);
albedo_snow = zeros(length(forz(:,1)), 1);
Mpot_snow = zeros(length(forz(:,1)), 1);
Mpot_ice = zeros(length(forz(:,1)), 1);
SWE = zeros(length(forz(:,1)), 1);
M_snow = zeros(length(forz(:,1)), 1);
EE= zeros(length(forz(:,1)), 1);
Qsnow = zeros(length(forz(:,1)), 1);
Qice = zeros(length(forz(:,1)), 1);
SWnet_snow=zeros(length(forz(:,1)), 1);
SWnet_ice=zeros(length(forz(:,1)), 1);
Pliq = zeros(length(forz(:,1)), 1);
Psnow = zeros(length(forz(:,1)), 1);
Em_snow = zeros(length(forz(:,1)), 1);
Em_ice = zeros(length(forz(:,1)), 1);
LW_net=zeros(length(forz(:,1)), 1);

% Iniciar variables de estado
area_gl(1) = Ai; % Se asigna el �rea inicial en km2
volumen(1) = cc * (Ai)^bb; % Se obtiene el volume inicial en km3
Sfree(1) = 0.5; % Se asigna el valor de Sfree inicial
albedo_snow(1) = 0.5; % albedo inicial

% parámetros
% días desde que cayó la ultima nevada
n = 0;

%% Análisis diario

Vice = 0;

ver = [];

conteo_agnos_SWE = 0;

SWE(1)=mean(mean_swe_year);

%% Modelo
for tt = 2:length(SWE)
    
    % 1. Se corrobora la precipitación
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t2m(tt) > Ts
        Pliq(tt) = pr(tt); % se asigna lo que cae como lluvia
    else
        Psnow(tt) = pr(tt); % se asigna lo que cae como nieve
        
        if pr(tt) > 5 % 5 [mm]/d
            n = 0; % se resetea el conteo de nieve
        else
            n = n+1;
        end
    end
    
    % 2. se calcula el albedo de la nieve
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        albedo_snow(tt) = albedo_snow_min + (albedo_snow_max - albedo_snow_min) * exp(-r * n); 
    
    if pr(tt) > 5 % 5 mm/d
        % se corrige como día con precipitación
        SWnet_snow(tt) = (1-albedo_snow(tt)) * sw(tt) * reduccion_forz_rain ; % [J/m2]
        SWnet_ice(tt) = (1-albedo_ice) * sw(tt) * reduccion_forz_rain ; % [J/m2]
    else
        SWnet_snow(tt) = (1-albedo_snow(tt)) * sw(tt) ; % [J/m2]
        SWnet_ice(tt) = (1-albedo_ice) * sw(tt); % [J/m2]
    end
    
    
     
     Em_ice(tt)=(SWnet_ice(tt)*sw_f+(t2m(tt)*C1)+C0)*(1000*3600*24/(Lf*rho_w));%+SWnet_snow/(Lf*rho_ice);
     Em_snow(tt)=(SWnet_snow(tt)*sw_f+(t2m(tt)*C1)+C0)*(1000*3600*24/(Lf*rho_w));%+SWnet_ice/(Lf*rho_snow)  
    
     Es_ice(tt)=(SWnet_ice(tt)*sw_f+(t2m(tt)*C1)+C0)*(1000*3600*24/(Ls*rho_w));%+SWnet_snow/(Lf*rho_ice);
     Es_snow(tt)=(SWnet_snow(tt)*sw_f+(t2m(tt)*C1)+C0)*(1000*3600*24/(Ls*rho_w));%+SWnet_ice/(Lf*rho_snow)
     
    %             clear SWnet_snow SWnet_ice
    
    % 4. Se calculan los derretimientos potenciales
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t2m(tt) >= To   % eliminate condition (stems from dergee day model? ???
        Mpot_snow(tt) = Em_snow(tt)*(1-Sfree(tt-1)); % [mm]
        Mpot_ice(tt) = rocoso * Sfree(tt-1) * Em_ice(tt); % [mm]
        
    else
        
        SUBLpot_snow(tt) = Es_snow(tt); % [mm]
        SUBLpot_ice(tt) = rocoso * Sfree(tt-1) * Es_ice(tt);
        
        Mpot_snow(tt) = 0 ;
        Mpot_ice(tt) = 0 ;
    end
    %
    if Mpot_snow(tt) < 0
        Mpot_snow(tt) = 0;    
    end
    if Mpot_ice(tt) < 0
        Mpot_ice(tt) = 0;  
    end
    
    if SUBLpot_snow(tt) < 0
        SUBLpot_snow(tt) = 0;    
    end
    if SUBLpot_ice(tt) < 0
        SUBLpot_ice(tt) = 0;  
    end
    
    
    % 5. Se obtiene el monto de SWE y el derretimiento
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SWE(tt) = SWE(tt-1) + Psnow(tt);
    M_snow(tt) = min(SWE(tt), Mpot_snow(tt));
    
    SUBL_snow= min(SWE(tt), SUBLpot_snow(tt));
    Sublimacion(tt)=SUBL_snow+SUBLpot_ice(tt);
    % 6. Obtener los caudales de la nieve y el hielo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qsnow(tt) = Qsnow(tt-1) * exp(-1/k_snow) + (Pliq(tt) + M_snow(tt)) * (1-exp(-1/k_snow)); % [mm]
    
    Qice(tt) = Qice(tt-1) * exp(-1/k_ice) + (Pliq(tt) * Sfree(tt) + Mpot_ice(tt)) * (1 - exp(-1/k_ice)); % [mm]
    
    % 7. Actualizar el monto de SWE sacando lo que se derriti� y Sfree
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SWE(tt) = SWE(tt) - M_snow(tt);
    
    if SWE(tt) < 0
        SWE(tt) = 0;
    end
    % Actualizar Sfree
    if Psnow(tt) > 0
        Sfree(tt) = 0;
    elseif SWE(tt) == 0
        Sfree(tt) = 1;
    elseif 1 - (SWE(tt) - M_snow(tt)) / SWE(tt) >= 0 & 1 - (SWE(tt) - M_snow(tt)) / SWE(tt) <= 1
        Sfree(tt) = 1 - (SWE(tt) - M_snow(tt)) / SWE(tt) ;
    else
        Sfree(tt) = 1;
    end
    % 8. Actualización del volumen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if tiempo_day2(tt,2) == actualizar_V(1) & tiempo_day2(tt,3) == actualizar_V(2)
        
        % Se incorpora el SWE si ha duforzo muchos años
        if SWE(tt) > 0
            conteo_agnos_SWE = conteo_agnos_SWE + 1;
        else
            conteo_agnos_SWE = 0;
        end
        
        if conteo_agnos_SWE == 1
            % Se vuelve a resetear el contador
            conteo_agnos_SWE = 0;
            
            volumen(tt) = volumen(tt-1) - Vice + calving + (1-Sfree(tt))*SWE(tt)*area_gl(tt-1)/10^6;
            % todo el SWE pasar a ser hielo
            SWE(tt) = 0;
        else
            volumen(tt) = volumen(tt-1) - Vice + calving;
            
        end
        
        ver = [ver Vice];
        
        
        % se escoge el "b" según el tamaño del glaciar
        %                     if area_gl(tt-1) > 30
        %                         bb = 1.3598;
        %                     else
        %                         bb = 1.27;
        %                     end
        
        area_gl(tt) = (volumen(tt) / cc ) ^ (1/bb);
        if volumen(tt) < 0
            volumen(tt) = 0;
            area_gl(tt) = 0;
        end
        Vice = 0;
    else
        area_gl(tt) = area_gl(tt-1);
        volumen(tt) = volumen(tt-1);
        Vice = Vice + ( area_gl(tt)*1000*1000 * Mpot_ice(tt)/1000 ) / 10^9;
    end
    
end

%Obtener el caudal efluente del glaciar
caudal_gl = 1000*Qice.*area_gl / (24*3600);
caudal_sn = 1000*Qsnow.*area_gl / (24*3600);
caudal_total = 1000*(Qice+Qsnow).*area_gl / (24*3600);
% Se obtiene un �ndice con años donde el SWE no es cero
filas = find(tiempo_day2(:,2)== 2 & tiempo_day2(:,3)==28);
SWE_index = [];
for jj = 2:length(filas)
    if min(SWE(filas(jj)-364:filas(jj))) == 0
        SWE_index = [SWE_index; 1];
    else
        SWE_index = [SWE_index; 0];
    end
end
C1_mod=pairs_zn(1,1);
C0_mod=pairs_zn(1,2);
flujos = [tiempo_day_output' area_gl volumen Sfree albedo_snow Mpot_snow Mpot_ice SWE M_snow Qsnow Qice caudal_gl Psnow caudal_sn caudal_total Sublimacion];
params = [Ts Ai cc bb a_ice sum(SWE_index) C1_mod C0_mod sw_f];

%% Modulo experimental para busqueda de predictores o variables influyentes en el derretimiento
tabla_cor=table();
tabla_cor.Msnow=M_snow;
tabla_cor.Mice=Mpot_ice;
tabla_cor.t2m=t2m;
tabla_cor.sw=sw;
tabla_cor.albedos=albedo_snow;
tabla_cor.SWE=SWE;
tabla_cor.swni=SWnet_ice;
tabla_cor.swns=SWnet_snow;
tabla_cor.pr=pr;
tabla_cor.prs=Psnow;
tabla_cor.area=area_gl;


end
