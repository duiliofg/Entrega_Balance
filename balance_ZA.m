function [params, flujos] = balance_ddi(rad, To, Ts, Ai, cc, bb, pairs, a_ice, rocoso, yr_mod, calving)
%% Parámetros
% rad: contiene la fecha, radiación de onda larga neta, onda corta
% incidente, precipitación diaria y temperatura [J/m2]
% To: Temperatura que separa la lluvia de la nieve [�C]
% r : tasa de decaimiento del Albedo [1/dias]
% Ai: �rea inicial del glaciar [km2]
% cc: Intercambio volumen glaciar
% bb: Intercambio volumen glaciar
% a_ice: albedo del glaciar
% Rf: factor que modifica la radiaci�n
tiempo_day2 = rad(1,1) : 1 : rad(end, 1);
tiempo_day_output=tiempo_day2;
tiempo_day2=datevec(tiempo_day2);
% fecha en la cual se realizará el cambio de
actualizar_V = [3, 15]; % [mes dia]
%%  Parámetros físicos
rho_w = 1000; % [kg/m3]
rho_ice=917; %[kg/m3]
rho_snow=400; %[kg/m3]
Lf = 334000; %[J/kg] Calor latente de fusión
B = 0.95; % Rango 0.8 a 1.1
Fice = 1; % Proporción del área superficial con hielo.
albedo_ice = a_ice; % albedo que tiene el hielo una vez que no hay nieve
k_ice=3.5; % factor que retrasa el derretimiento del hielo
k_snow=12.5; % factor que retrasa el derretimeinto de la nieve
albedo_snow_min = 0.4;
albedo_snow_max = 0.85;
reduccion_rad_rain = 0.3; % reducción de la radiación en días con lluvia DeWalle&Rango

% Caso de calving
if isnan(calving)==true
    calving=0;
else 
    calving=calving;
end 

% Crear vectores con variables de estado
area_gl = zeros(length(rad(:,1)), 1);
volumen = zeros(length(rad(:,1)), 1);
Sfree = zeros(length(rad(:,1)), 1);
albedo_snow = zeros(length(rad(:,1)), 1);
Mpot_snow = zeros(length(rad(:,1)), 1);
Mpot_ice = zeros(length(rad(:,1)), 1);
SWE = zeros(length(rad(:,1)), 1);
M_snow = zeros(length(rad(:,1)), 1);
Qsnow = zeros(length(rad(:,1)), 1);
Qice = zeros(length(rad(:,1)), 1);
Pliq = zeros(length(rad(:,1)), 1);
Psnow = zeros(length(rad(:,1)), 1);
Em_snow = zeros(length(rad(:,1)), 1);
Em_ice = zeros(length(rad(:,1)), 1);

% Iniciar variables de estado
area_gl(1) = Ai; % Se asigna el área inicial en km2
volumen(1) = cc * (Ai)^bb; % Se obtiene el volume inicial en km3
Sfree(1) = 0.5; % Se asigna el valor de Sfree inicial
albedo_snow(1) = 0.5; % albedo inicial



% días desde que cay� la ultima nevada
n = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análisis diario
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vice = 0;

ver = [];

conteo_agnos_SWE = 0;

%SWE calculo inicial

for xx = 2:3653

    if rad(xx,3) > Ts
        Pliq(xx) = rad(xx,2); % se asigna lo que cae como lluvia
    else
        Psnow(xx) = rad(xx,2); % se asigna lo que cae como nieve

        if rad(xx,2) > 5 % 5 [mm]/d
            n = 0; % se resetea el conteo de nieve
        else
            n = n+1;
        end
    end


    % 3.1 Degree-day
    if rad(xx,3)>0
        H(xx)=1;
    else
        H(xx)=0;
    end
    D(xx)=rad(xx,3)*H(xx);
    Em_ice(xx)=D(xx)*pairs(1,1);%+SWnet_snow/(Lf*rho_ice);
    Em_snow(xx)=D(xx)*pairs(1,2);%+SWnet_ice/(Lf*rho_snow);

    %             clear SWnet_snow SWnet_ice

    % 4. Se calculan los derretimientos potenciales
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rad(xx,3) > To
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


    % 5. Se obtiene el monto de SWE y el derretimiento
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SWE(xx) = SWE(xx-1) + Psnow(xx);
    M_snow(xx) = min(SWE(xx), Mpot_snow(xx));

    % 6. Obtener los caudales de la nieve y el hielo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qsnow(xx) = Qsnow(xx-1) * exp(-1/k_snow) + (Pliq(xx) + M_snow(xx)) * (1-exp(-1/k_snow)); % [mm]

    Qice(xx) = Qice(xx-1) * exp(-1/k_ice) + (Pliq(xx) * Sfree(xx) + Mpot_ice(xx)) * (1 - exp(-1/k_ice)); % [mm]

    % 7. Actualizar el monto de SWE sacando lo que se derriti� y Sfree
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    % 8. Actualización del volumen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if tiempo_day2(xx,2) == actualizar_V(1) & tiempo_day2(xx,3) == actualizar_V(2)

        % Se incorpora el SWE si ha durado muchos años
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

area_gl = zeros(length(rad(:,1)), 1);
volumen = zeros(length(rad(:,1)), 1);
Sfree = zeros(length(rad(:,1)), 1);
albedo_snow = zeros(length(rad(:,1)), 1);

Mpot_snow = zeros(length(rad(:,1)), 1);
Mpot_ice = zeros(length(rad(:,1)), 1);
SWE = zeros(length(rad(:,1)), 1);
M_snow = zeros(length(rad(:,1)), 1);

Qsnow = zeros(length(rad(:,1)), 1);
Qice = zeros(length(rad(:,1)), 1);

Pliq = zeros(length(rad(:,1)), 1);
Psnow = zeros(length(rad(:,1)), 1);

Em_snow = zeros(length(rad(:,1)), 1);
Em_ice = zeros(length(rad(:,1)), 1);

% Iniciar variables de estado
area_gl(1) = Ai; % Se asigna el �rea inicial en km2
volumen(1) = cc * (Ai)^bb; % Se obtiene el volume inicial en km3
Sfree(1) = 0.5; % Se asigna el valor de Sfree inicial
albedo_snow(1) = 0.5; % albedo inicial

% parámetros
% días desde que cayó la ultima nevada
n = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Análisis diario
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vice = 0;

ver = [];

conteo_agnos_SWE = 0;

SWE(1)=mean(mean_swe_year);

%% Modelo
for tt = 2:length(SWE)

    % 1. Se corrobora la precipitación
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rad(tt,3) > Ts
        Pliq(tt) = rad(tt,2); % se asigna lo que cae como lluvia
    else
        Psnow(tt) = rad(tt,2); % se asigna lo que cae como nieve

        if rad(tt,2) > 5 % 5 [mm]/d
            n = 0; % se resetea el conteo de nieve
        else
            n = n+1;
        end
    end


    % 3.1 Degree-day (Autor: Duilio Fonseca)
    if rad(tt,3)>0
        H=1;
    else
        H=0;
    end
    D=rad(tt,3)*H;
    Em_ice(tt)=D*pairs(1,1);
    Em_snow(tt)=D*pairs(1,2);

    %             clear SWnet_snow SWnet_ice

    % 4. Se calculan los derretimientos potenciales
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rad(tt,3) > To
        Mpot_snow(tt) = Em_snow(tt); % [mm]
        Mpot_ice(tt) = rocoso * Sfree(tt-1) * Em_ice(tt); % [mm]
    else
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




    % 5. Se obtiene el monto de SWE y el derretimiento
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SWE(tt) = SWE(tt-1) + Psnow(tt);
    M_snow(tt) = min(SWE(tt), Mpot_snow(tt));

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

        % Se incorpora el SWE si ha durado muchos años
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
DD_I_mod=pairs(1,1);
DD_S_mod=pairs(1,2);
flujos = [tiempo_day_output' area_gl volumen Sfree albedo_snow Mpot_snow Mpot_ice SWE M_snow Qsnow Qice caudal_gl Psnow caudal_sn caudal_total];
params = [Ts Ai cc bb a_ice sum(SWE_index) DD_I_mod DD_S_mod];

end
