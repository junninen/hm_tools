
e=event_table_load_all();
Iev=find(e.h.class1+e.h.class2);
% st=[];
% nd=[];
% mxN=[];
for i=137:length(Iev)
    
    hm_plot(hm_load(e.h.daten(Iev(i)),[],'dmps')),caxis([1 inf])
    gout=ginput(3);
    
    st(i,:)=[e.h.daten(Iev(i)),gout(1,:)];
    nd(i,:)=[e.h.daten(Iev(i)),gout(2,:)];
    mxN(i,:)=[e.h.daten(Iev(i)),gout(3,:)];
    
end