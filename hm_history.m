function hmD=hm_history(hmD,command)
%
%Write command history that have modified the hm-structure
%
% hmD=hm_history(hmD,command)
%
% Modifies hmD.meta.history
% hmD.meta.history = 2xn cell, time of action, command used
%

%Heikki Junninen
%Nov 2007

hmD.meta.history=[hmD.meta.history;{datestr(now,30),command}];