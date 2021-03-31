function f_creer_newvar2(nc,namevar,type,dimid,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cascade Exploitation :
% -----------------------
%
%  Fonction permettant de creer une nouvelle 
%   Variable dans un fichier netcdf deja ouvert.
%   Fonction a utiliser pour creer plusieurs variables dans un meme fichier NetCDF existant.
%
%   f_creer_newvar2(filenc,namevar,type,dimvar,varargin)
%
% En entree :
% -------------
%   nc : identifiant du fichier NetCDF ou la variable doit etre creee.
%   namevar : nom de la variable a creer
%   type : Type de la variable
%   dimid :  Liste des identifiants des dimensions de la variable. 
%   varargin : liste(nom,valeur) des attributs de la variable a creer.  ex: 'long_name','Vitessse en ms'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%      C. Kermabon - P Le Bot Avril 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition/Creation de la variable.
id_var=netcdf.defVar(nc,namevar,type,dimid);
% Et de ses attributs.
for iatt=1:2:size(varargin,2)
 netcdf.putAtt(nc,id_var,varargin{iatt},varargin{iatt+1});
end
