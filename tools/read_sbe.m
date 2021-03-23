function [data,tmp1_sensor_scale,tmp2_sensor_scale,prs_sensor_type,oxy1_sensor_type,oxy2_sensor_type,turb_sensor_type,trans_sensor_type,fluo_sensor_type,cpar_sensor_type] = read_sbe(fname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fonction matlab permettant de lire les fichiers ASCII *.cnv
% issus de l'acquisition SEABIRD.
% Les fichiers *CNV contiennent une entete suivie des donnees.
%
% En entree : fname = nom du fichier ASCII Seabird a lire
% En sortie : data = Structure contenant les donnees CTD du fichier ainsi que
%             les donnees de position et de jour.
%             tmp1_sensor_scale: Echelle de mesure du premier capteur de temperature  (ex: ITS-90).
%             tmp2_sensor_scale: Echelle de mesure du second capteur de temperature  (ex: ITS-90).
%             prs_sensor_type: Type du capteur de pression (ex: Digiquartz)
%             oxy1_sensor_type : Type du capteur premier capteur d'oxygene (ex: SBE43)
%             oxy2_sensor_type : Type de capteur du second capteur d
%             oxygene (ex: SBE43)
%
% C. Kermabon ; 15/03/2011
%
% P. Le Bot : 26/09/2011 - Ajout.
%   Recupï¿½ration des infos type de capteur Prs, Oxy1, Oxy2.
%   Recuperation de l'echelle de mesure du capteur de temperature (ex. ITS-90).
% Modifs Pierre Rousselot - US191 IMAGO, 03/2017
%       Ajout fluorimetre et transmissiometre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(fname,'r');
if fid<0
    h_dlg = warndlg(sprintf('Attention\nFichier ''%s'' inexistant\n',fname));
    waitfor(h_dlg);
    return;
end
%
% Lecture de l'entete.
%
tline = fgets(fid);
data.entete = [];
liste_param = '';
num_param = [];
unit_param = '';

% On initialise les infos capteur
tmp1_sensor_scale = 'None';
tmp2_sensor_scale = 'None';
prs_sensor_type   = 'None';
oxy1_sensor_type  = 'None';
oxy2_sensor_type  = 'None';
turb_sensor_type  = 'None';
trans_sensor_type = 'None';
fluo_sensor_type  = 'None';
cpar_sensor_type  = 'None';

while isempty(strfind(tline,'*END'))
    data.entete = char(data.entete,tline);
    % Recuperation du nombre de parametre du fichier.
    if (~isempty(strfind(tline,'nquan')))
        nparam = textscan(tline,'%*s %*s %*s %d','delimiter',' ');
        nparam = nparam{1}(1);
        if exist('param_PMEL', 'var')
            nparam = 21;
        end
    end
    % Recuperation des colonnes associees aux divers parametres CTD
    % ainsi qu'a la latitude, la longitude, le jour julien et l'annee de la CTD.
    if (~isempty(strfind(tline,'name')))
        if (~isempty(strfind(tline,'prDM')))
            liste_param = char(liste_param,'prs');
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %s %s');
            prs_sensor_type=char(strcat(sensor_type{1},sensor_type{2}));
            ind_pression = 1;
        elseif (~isempty(strfind(tline,'t0')))
            liste_param = char(liste_param,'temp1');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_scale= textscan(tline,'%*s %*s %*d %*s %*s %*s [%s, %*s');
            tmp1_sensor_scale=char(sensor_scale{1});
            if (strcmp(tmp1_sensor_scale(end),','))
                tmp1_sensor_scale=tmp1_sensor_scale(1:end-1);
            end
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
            ind_temp1 = 1;
        elseif (~isempty(strfind(tline,'t1')))
            liste_param = char(liste_param,'temp2');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_scale= textscan(tline,'%*s %*s %*d %*s %*s %*s %*d [%s, %*s');
            tmp2_sensor_scale=char(sensor_scale{1});
            if (strcmp(tmp2_sensor_scale(end),','))
                tmp2_sensor_scale=tmp2_sensor_scale(1:end-1);
            end
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
            ind_temp2 = 1;
        elseif (~isempty(strfind(tline,'c0')))
            liste_param = char(liste_param,'cond1');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param, ' ');
            end
            ind_cond1 = 1;
        elseif(~isempty(strfind(tline,'c1')))
            liste_param = char(liste_param,'cond2');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
            ind_cond2 = 1;
        elseif(~isempty(strfind(tline,'sbeox0V')))
            liste_param = char(liste_param,'oxy1');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s,');
            oxy1_sensor_type=char(strcat(sensor_type{1},sensor_type{2}));
            if (strcmp(oxy1_sensor_type(end),','))
                oxy_sensor_type=oxy1_sensor_type(1:end-1);
            end
            unit_param = char(unit_param,'Volt');
            ind_oxy1 = 1;
        elseif (~isempty(strfind(tline,'sbeox1V')))
            liste_param = char(liste_param,'oxy2');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s,');
            oxy2_sensor_type=char(strcat(sensor_type{1},sensor_type{2}));
            if (strcmp(oxy2_sensor_type(end),','))
                oxy2_sensor_type=oxy2_sensor_type(1:end-1);
            end
            unit_param = char(unit_param,'[Volt]');
            ind_oxy2 = 1;
        elseif(~isempty(strfind(tline,'sbeox0M')))
            liste_param = char(liste_param,'oxyl1');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s,');
            unit_param = char(unit_param,'[Volt]');
            ind_oxy1 = 1;
        elseif (~isempty(strfind(tline,'sbeox1M')))
            liste_param = char(liste_param,'oxyl2');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s,');
            oxy2_sensor_type=char(strcat(sensor_type{1},sensor_type{2}));
            if (strcmp(oxy2_sensor_type(end),','))
                oxy2_sensor_type=oxy2_sensor_type(1:end-1);
            end
            unit_param = char(unit_param,'Volt');
            ind_oxy2 = 1;
        elseif (~isempty(strfind(tline,'ptemp')))
            liste_param = char(liste_param,'ptemp');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'latitude')))
            liste_param = char(liste_param,'latitude');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_latitude = 1;
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'longitude')))
            liste_param = char(liste_param,'longitude');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'Julian Days')))
            liste_param = char(liste_param,'jourjul');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_jour = 1;
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'flECO-AFL')))
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s,');
            if strcmp(fluo_sensor_type,'None')
                liste_param = char(liste_param,'fluo');
                fluo_sensor_type={char(strcat(sensor_type{1},sensor_type{2}))};
            else
                liste_param = char(liste_param,'fluo2');
                fluo_sensor_type=[{fluo_sensor_type} {char(strcat(sensor_type{1},sensor_type{2}))}];
            end                       
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'flC')))
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s,');
            if strcmp(fluo_sensor_type,'None')
                liste_param = char(liste_param,'fluo');
                fluo_sensor_type={char(strcat(sensor_type{1},sensor_type{2}))};
            else
                liste_param = char(liste_param,'fluo2');
                fluo_sensor_type=[{fluo_sensor_type} {char(strcat(sensor_type{1},sensor_type{2}))}];
            end            
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'CStarTr0')))         
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s %s,');
            if strcmp(trans_sensor_type,'None')
                liste_param = char(liste_param,'trans');
                trans_sensor_type={char(strcat(sensor_type{1},sensor_type{2},sensor_type{3}))};
            else
                liste_param = char(liste_param,'trans2');
                trans_sensor_type=[{trans_sensor_type} {char(strcat(sensor_type{1},sensor_type{2},sensor_type{3}))}];                
            end
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'CStarTr1')))
            liste_param = char(liste_param,'trans2');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %*s %s %s %s,');
            sens_type=char(strcat(sensor_type{1},sensor_type{2},sensor_type{3}));
            trans_sensor_type=[trans_sensor_type;sens_type(1:end-1)];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'cpar')))
            liste_param = char(liste_param,'cpar');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %s %s,');
            cpar_sensor_type=char(strcat(sensor_type{1},sensor_type{2}));
            if (strcmp(cpar_sensor_type(end),','))
                cpar_sensor_type=cpar_sensor_type(1:end-1);
            end
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'sigma-é00')))
            liste_param = char(liste_param,'dens1');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'sigma-é11')))
            liste_param = char(liste_param,'dens2');
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'turbWETntu0')))
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %s %s %s %s,');
            if strcmp(turb_sensor_type,'None')
                liste_param = char(liste_param,'turb');
                turb_sensor_type=[char(strcat(sensor_type{1},sensor_type{2},sensor_type{3})) '  '];
            else
                liste_param = char(liste_param,'turb2');
                turb_sensor_type=[turb_sensor_type [char(strcat(sensor_type{1},sensor_type{2},sensor_type{3})) '  ']];
            end
            
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        elseif (~isempty(strfind(tline,'turbWETbb0')))
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            sensor_type=textscan(tline,'%*s %*s %*d %*s %*s %*s %s %s %s %s,');
            if strcmp(turb_sensor_type,'None')
                liste_param = char(liste_param,'turb');
                turb_sensor_type=char(strcat(sensor_type{1},sensor_type{2},sensor_type{3},sensor_type{4}));
            else
                liste_param = char(liste_param,'turb2');
                turb_sensor_type=[turb_sensor_type;char(strcat(sensor_type{1},sensor_type{2},sensor_type{3},sensor_type{4}))];
            end                       
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param = char(unit_param,' ');
            end
        else
            ind_egal = strfind(tline,'=');
            ind_point = strfind(tline,':');
            nom_bid = tline(ind_egal+2:ind_point-1);
            ind_slash = strfind(nom_bid,'/');
            nom_bid(ind_slash) ='_';
            liste_param = char(liste_param,nom_bid);
            num_param = [num_param textscan(tline,'%*s %*s %d %*d','delimiter',' ')];
            ind_crochet = strfind(tline,'[');
            ind_crochet2 = strfind(tline,']');
            if ~isempty(ind_crochet)
                unit_param = char(unit_param,tline([ind_crochet:ind_crochet2]));
            else
                unit_param= char(unit_param,' ');
            end
        end
    end
    if (~isempty(strfind(tline,'NMEA UTC')))
        data.an_ctd = sscanf(tline,'%*c %*s %*s %*s  = %*s %*d %d');
    end
    if (~isempty(strfind(tline,'bad_flag')))
        data.badval = sscanf(tline,'%*c %*s = %f');
    end
    tline = fgets(fid);
end
data.entete = char(data.entete,tline);
data.entete = data.entete(2:end,:);
%
% Lecture des donnees.
%
chaine_format = '';
for i_param=1:nparam
    chaine_format = [chaine_format ' %f'];
end
res = fscanf(fid,chaine_format,[nparam Inf]);
fclose(fid);

for i_param = 1:length(num_param)
    chaine = sprintf('data.%s=res(%d,:);',strtrim(liste_param(i_param+1,:)),num_param{i_param}+1);
    eval(chaine);
end
data.liste_param = liste_param(2:end,:);
data.unit_param = unit_param(2:end,:);


