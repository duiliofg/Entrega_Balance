function [min_dhdt,best_dhdt,diff_dhdt,dhdt_val,diff_dhdt_noabs]=validacion_dhdt_hugonnet(tiempo_day,dhdt_val,flujos_acum,yr_mod)
%dhdt_val: variación de volumen por año según producto de Hugonnet et al. (2021)
%Resultados acumulados para glaciar.
%% Revisión de meteodología de selección. 
% Parametros temporales

periodo_de_validacion=(datenum(yr_mod,03,16,00,00,00):datenum(0000,00,01,00,00,00):datenum(2019,12,31,00,00,00))';
index_final=find(tiempo_day==datenum(2019,12,31));
vector_tiempo=datevec(periodo_de_validacion);
yr=vector_tiempo(:,1);
mth=vector_tiempo(:,2);
day=vector_tiempo(:,3);
yr_frac=(periodo_de_validacion(end)-periodo_de_validacion(1))/365.25;
% Carga de tablas para cálculo de balance de masa
MB=table();
MB.yr=yr;
MB.mth=mth;
MB.day=day;
dhdt_val=(917/1000)*dhdt_val; % conversión de [m.ice.eq] a [m.w.eq]

for tt=1:size(flujos_acum,3)
muestra=tt; % indice representa una modelacion bajo una parametrización determinada (ver pair_zn)
MB.area=flujos_acum(1:index_final,2,tt);
MB.volumen=flujos_acum(1:index_final,3,tt);
dhdt_modelado(tt,1)=(MB.volumen(end)-MB.volumen(1))*1000/MB.area(1)/yr_frac;

% liberando variables
MB.area=[];
MB.volumen=[];

clear MB_year

end
diff_dhdt_r=dhdt_val-dhdt_modelado;
diff_dhdt=abs(dhdt_val-dhdt_modelado);
min_dhdt=find(diff_dhdt==min(diff_dhdt));

if length(min_dhdt)>1
    min_dhdt=min_dhdt(1);
else
    %Nothing
end 

diff_dhdt_noabs=diff_dhdt_r(min_dhdt);

best_dhdt=dhdt_modelado(min_dhdt);
diff_dhdt=diff_dhdt(min_dhdt);

end
