function i_str=hm_readInstrumentParam(fName)
%
%reads parameters from HMTools parameters file
%
%File structure must be
%<instrument>
%  <name>SomeName</name>
%  <path>SomePath</path>
%  <fnm>sprintf(SomeModel)</fnm>
%  <alt_fnm>sprintf(SomeModel)</alt_fnm> %optional
%</instrument>
%<instrument>
%  <name>SomeOtherName</name>
%  <path>SomeOtherPath</path>
%  <fnm>sprintf(SomeOtherModel)</fnm>
%</instrument>
%
%used by hm_load

%Heikki Junninen 20.02.2007

fid=fopen(fName,'r');

i_col=0; %switch for instrument collector, 1=collect 0=dont
name_col=0; %name collector switch
path_col=0; %path collector switch
fnm_col=0; %fnm collector switch
alt_fnm_col=0; %fnm collector switch
nrInstr=0; %number of instruments counter

while 1
    tline=fgetl(fid);
    if tline==-1, break,end
    tline = deblank(tline); %remove tailing blanks
    tline = fliplr(deblank(fliplr(tline))); %remove front blanks
    
    if ~isempty(tline) && tline(1)~='%',
   
        %switch instrument collector on or off
        if strfind(tline,'<instrument>'),
            i_col=1;
            nrInstr=nrInstr+1;
        end
        
        if strfind(tline,'</instrument>'),
            i_col=0;
        end
        
        %collect name
        %value on different row than tags
        if i_col & name_col,
            i_str(nrInstr).name=tline;
            name_col=0;
        end
        if i_col & strfind(tline,'<name>'),
            name_col=1;
        end
        
         %tags on the same row with value       
        if i_col & strfind(tline,'<name>') & strfind(tline,'</name>'),
            i_str(nrInstr).name=tline(strfind(tline,'<name>')+6:strfind(tline,'</name>')-1);
            name_col=0;
        end
 
        %collect path
        %value on different row than tags
        if i_col & path_col,
            i_str(nrInstr).path=tline;
            path_col=0;
        end
        if i_col & strfind(tline,'<path>'),
            path_col=1;
        end
        
         %tags on the same row with value       
        if i_col & strfind(tline,'<path>') & strfind(tline,'</path>'),
            i_str(nrInstr).path=tline(strfind(tline,'<path>')+6:strfind(tline,'</path>')-1);
            path_col=0;
        end

        %collect fnm
        %value on different row than tags
        if i_col & fnm_col,
            i_str(nrInstr).fnm=tline;
            fnm_col=0;
        end
        if i_col & strfind(tline,'<fnm>'),
            fnm_col=1;
        end
        
         %tags on the same row with value       
        if i_col & strfind(tline,'<fnm>') & strfind(tline,'</fnm>'),
            i_str(nrInstr).fnm=tline(strfind(tline,'<fnm>')+5:strfind(tline,'</fnm>')-1);
            fnm_col=0;
        end
               
        %collect fnm
        %value on different row than tags
        if i_col & alt_fnm_col,
            i_str(nrInstr).alt_fnm=tline;
            alt_fnm_col=0;
        end
        if i_col & strfind(tline,'<alt_fnm>'),
            alt_fnm_col=1;
        end
        
         %tags on the same row with value       
        if i_col & strfind(tline,'<alt_fnm>') & strfind(tline,'</alt_fnm>'),
            i_str(nrInstr).alt_fnm=tline(strfind(tline,'<alt_fnm>')+9:strfind(tline,'</alt_fnm>')-1);
            alt_fnm_col=0;
        end

    end
end
fclose(fid);