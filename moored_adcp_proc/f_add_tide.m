function f_add_tide(filenc)
% ==========================================================
%
% Cascade Exploitation :
% ------------------------
% Fonction permettant de rajouter les informations de maree 
% dans le fichier NetCDF en cours de traitement.
%
% En entree :
% ------------
% - filenc : nom du fichier NetCDF a traiter.
%
% ==========================================================
%
% 2009 : P.Le Bot - C. Kermabon.
%
% ==========================================================
%global_CE;        
%hbar = waitbar(0,'Tide : 11 harmonics used'); % Barre d'evolution du travail en cours.
%
% Récuperation des dates et positions des données ADCP.
%
%nc=netcdf.open(filenc,'NOWRITE');
%tpsID=netcdf.inqVarID(nc,'JULD');
%tps=netcdf.getVar(nc,tpsID);
%drefID=netcdf.inqVarID(nc,'REFERENCE_DATE_TIME');
%dref=netcdf.getVar(nc,drefID);
%if (size(dref,1)~=1)
%         dref=dref';
%end
%netcdf.close(nc);
time= data.time;
lat=mooring.lat*ones(length(time),1);%f_autonan(filenc,'LATITUDE');
lon=mooring.lon*ones(length(time),1);%f_autonan(filenc,'LONGITUDE');
u_interp = data.u_interp;
v_interp = data.v_interp;
%tps= tps+ (jul_0h(str2double(dref(1:4)),str2double(dref(5:6)),str2double(dref(7:8))));
%res=greg_0h(tps);
%tps_tide = datenum(res); % Pour l'appel à tide_pred, les dates doivent etre fournies sous forme datenum.
for k=1:length(time)
    [yr, mn, dy, hr]= gregorian(time(k));
    tps_tide(k)=datenum(yr, mn, dy, hr, 00, 00);
end
%
% Récuperation du modele TPX utilisé fourni par Egbert and Erofeeva.
% Puis appel aux subroutines de calcul associées (cf.  http://www.esr.org/polar_tide_models/Model_TPXO71.html)
model_maree = 'C:\Workspace_Matlab\ADCP_mooring\moored_adcp_proc\tide\model\Model_tpxo7.1';%%fullfile(getenv('REP_MAREE'),filesep,getenv('MODEL_MAREE'));
% Calcul des vitesses de maree (cm/s).
nb_valeur = length(tps_tide);%length(lat);
U_maree = nan*ones(nb_valeur,1);
V_maree = nan*ones(nb_valeur,1);
U_pred = nan*ones(nb_valeur,1);
V_pred = nan*ones(nb_valeur,1);
isok = find(isfinite(tps_tide));%find(isfinite(lat));
U_maree(isok) = tide_pred(model_maree,tps_tide(isok),lat(isok),lon(isok),'u',[1:11]); % Calcul avec les 11 harmoniques principales.
%waitbar(0.25,hbar);
V_maree(isok) = tide_pred(model_maree,tps_tide(isok),lat(isok),lon(isok),'v',[1:11]);
%waitbar(0.50,hbar);
% Calcul des transport de maree (m^2/s).
U_pred(isok) = tide_pred(model_maree,tps_tide(isok),lat(isok),lon(isok),'U',[1:11]);
%waitbar(0.75,hbar);
V_pred(isok) = tide_pred(model_maree,tps_tide(isok),lat(isok),lon(isok),'V',[1:11]);
%waitbar(1,hbar);

% Ajout des variables associees a la maree dans le fichier.
% ---------------------------------------------------------
% Si les variables associees a la maree n'existent pas, on les cree.
% % 
% res_var=f_test_vars(filenc,'U_TIDE');
% if res_var==0
%        nc = netcdf.open(filenc,'NC_WRITE');
%        netcdf.reDef(nc);
%        dim0 = netcdf.inqDimID(nc,'N_DATE_TIME');
%        dim1 = netcdf.inqDimID(nc,'N_LEVEL');
%        f_creer_newvar2(nc,'U_TIDE','NC_FLOAT',dim0,'units','meter per second','type_tide',getenv('MODEL_MAREE'),'long_name','Eastward tide Velocity','_FillValue',single(-999999),'valid_min',-20,'valid_max',20);
%        f_creer_newvar2(nc,'V_TIDE','NC_FLOAT',dim0,'units','meter per second','type_tide',getenv('MODEL_MAREE'),'long_name','Northward tide Velocity','_FillValue',single(-999999),'valid_min',-20,'valid_max',20);
%        f_creer_newvar2(nc,'UVEL_ADCP_CORTIDE','NC_FLOAT',[dim1 dim0],'units','meter per second','long_name','Eastward absolute velocity corrected for tide','_FillValue',single(-999999),'valid_min',-20,'valid_max',20); 
%        f_creer_newvar2(nc,'VVEL_ADCP_CORTIDE','NC_FLOAT',[dim1 dim0],'units','meter per second','long_name','Northward absolute velocity corrected for tide','_FillValue',single(-999999),'valid_min',-20,'valid_max',20); 
%        netcdf.endDef(nc);
%        netcdf.close(nc);
% end
% res_var=f_test_vars(filenc,'TU_TIDE');
% if res_var==0
%        nc = netcdf.open(filenc,'NC_WRITE');
%        netcdf.reDef(nc);
%        dim0 = netcdf.inqDimID(nc,'N_DATE_TIME');
%        dim1 = netcdf.inqDimID(nc,'N_LEVEL');
%        f_creer_newvar2(nc,'TU_TIDE','NC_FLOAT',dim0,'units','meter^2 per second','long_name','Eastward tide Transport','_FillValue',single(-999999),'valid_min',-20,'valid_max',20);
%        f_creer_newvar2(nc,'TV_TIDE','NC_FLOAT',dim0,'units','meter^2 per second','long_name','Northward tide Transport','_FillValue',single(-999999),'valid_min',-20,'valid_max',20);  
%        netcdf.endDef(nc);
%        netcdf.close(nc);
% end
% 
% 
% 
% % Ecriture des valeurs de vitesses de maree.
% nc = netcdf.open(filenc,'NC_WRITE');
% 
% % Mise a jour de l'attribut du nom de maree utilise (type_tide)
% netcdf.reDef(nc);
% netcdf.putAtt(nc,netcdf.inqVarID(nc,'U_TIDE'),'type_tide',getenv('MODEL_MAREE'));
% netcdf.putAtt(nc,netcdf.inqVarID(nc,'V_TIDE'),'type_tide',getenv('MODEL_MAREE'));
% netcdf.endDef(nc);
% 
% idv = netcdf.inqVarID(nc,'V_TIDE');
% idu = netcdf.inqVarID(nc,'U_TIDE');                      
U_maree=U_maree/100;  % Passage de cm/s en m/s
V_maree=V_maree/100;
indnan= isnan(U_maree);  
figure;
hist(U_maree,100); title('U maree (m/s)');
figure;
hist(V_maree,100); title('V maree (m/s)');

% U_maree(indnan)=f_get_fillvalue(filenc,'U_TIDE');
% clear indnan;
% indnan= isnan(V_maree);
% V_maree(indnan)=f_get_fillvalue(filenc,'V_TIDE'); 
% netcdf.putVar(nc,idu,U_maree);
% netcdf.putVar(nc,idv,V_maree);
% % Idem pour le Transport maree.   
% idv = netcdf.inqVarID(nc,'TV_TIDE');
% idu = netcdf.inqVarID(nc,'TU_TIDE');          
% indnan= isnan(U_pred);
% U_pred(indnan)=f_get_fillvalue(filenc,'TU_TIDE');
% clear indnan;
% indnan= isnan(V_pred);
% V_pred(indnan)=f_get_fillvalue(filenc,'TV_TIDE');
% clear indnan;
% netcdf.putVar(nc,idu,U_pred(1:nb_valeur));
% netcdf.putVar(nc,idv,V_pred(1:nb_valeur));
% ====================================
%   Calcul des vitesses corrigees:
% ====================================
% Lecture vitesses U.
% -----------------------
figure;
subplot(1,2,1);
hist(u_interp(:),100);
Nbins=size(u_interp,2);
%
% On corrige les vitesses uniquement lorsque la maree a pu etre
% calculee (ie different de -999999 et Inf). Sinon, on garde la vitesse non corrigee
% de la maree.
%
isok = find( (isfinite(U_maree)));% & (U_maree~=f_get_fillvalue(filenc,'U_TIDE')));
u_interp_cor(isok,:) = u_interp(isok,:)-U_maree(isok)*ones(1,Nbins); 

subplot(1,2,2);
hist(u_interp_cor(:),100);

% Ecriture de V corrigee.
%--------------------------
figure;
subplot(1,2,1);
hist(v_interp(:),100);
Nbins=size(v_interp,2);
%
% On corrige les vitesses uniquement lorsque la maree a pu etre
% calculee (ie different de -999999 et Inf). Sinon, on garde la vitesse non corrigee
% de la maree.
%
isok = find( (isfinite(V_maree)));% & (U_maree~=f_get_fillvalue(filenc,'U_TIDE')));
v_interp_cor(isok,:) = v_interp(isok,:)-V_maree(isok)*ones(1,Nbins); 

subplot(1,2,2);
hist(v_interp_cor(:),100);

%netcdf.putVar(nc,idu,Vit_cor(:,:)');
%
% Ecriture de V corrigee.
%
% Vit_cor=f_autonan(filenc,'VVEL_ADCP')';
% Nbins=size(Vit_cor,2);
%
% On corrige les vitesses uniquement lorsque la maree a pu etre
% calculee (ie different de -999999 et Inf). Sinon, on garde la vitesse non corrigee
% de la maree.
%
% isok = find( (isfinite(V_maree)) & (V_maree~=f_get_fillvalue(filenc,'V_TIDE')));
% Vit_cor(isok,:) = Vit_cor(isok,:)-V_maree(isok)*ones(1,Nbins); 
% %idv = netcdf.inqVarID(nc,'VVEL_ADCP_CORTIDE');
% clear indnan;
% indnan= isnan(Vit_cor);
% Vit_cor(indnan)=f_get_fillvalue(filenc,'VVEL_ADCP_CORTIDE');
%netcdf.putVar(nc,idv,Vit_cor(:,:)');

% Mise a jour de l'attribut DATE_UPDATE
%  netcdf.putAtt(nc,netcdf.getConstant('NC_GLOBAL'),'DATE_UPDATE',date);
% 
% netcdf.close(nc);
% delete(hbar);




