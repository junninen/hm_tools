function [str_out,pos]=hm_replace_char(str,txt1,txt2)
%
% replace txt1 in str with txt2
% str_out=hm_replace_char(str,txt1,txt2)
%
%

%Heikki Junninen 18.05.2007


Itxt1=strfind(str,txt1); %position of txt1
pos=Itxt1;
Ltxt1=length(txt1); %length of txt1

while ~isempty(Itxt1)
    str1=str(1:max(Itxt1(1)-1,1));
    str2=str(Itxt1(1)+Ltxt1:end);
    if Itxt1(1)==1
        str=[txt2,str2];
    else
        str=[str1,txt2,str2];
    end
    Itxt1=strfind(str,txt1); %find again since length of str might have changed if many txt1 found
end

str_out=str;
