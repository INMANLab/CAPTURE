function RegionTable = ExtractRegionNames(ChanLabels)

chNum = size(ChanLabels,2);
pNum = size(ChanLabels,1);
%Finding regions
RegionTable = [];
for pIdx=1:pNum
    for chIdx=1:chNum %chan
        chanlabel = ChanLabels{pIdx,chIdx};
        b1 = contains(chanlabel,("Amygdala"|"Anterior Hippocampus"));
        b2 = contains(chanlabel,("Temporal"));
        b3 = contains(chanlabel,("Entorhinal"|"Perirhinal"));
        if b1
            llb = 1; %anterior
            llbname = 'AntHipp';
        else %posterior
            if b2 %temporal
                llb = 2;
                llbname = 'LatTemp';
            elseif b3 %entorhinal
                llb = 3;
                llbname = 'Ent+Peri';
            else
                llb = 4;
                llbname = 'PostHipp+Para';
            end
        end
        RegionTable = vertcat(RegionTable,table(pIdx,chIdx,{chanlabel},llb,{llbname},'VariableNames',{'Patient','Chan','ChanLabel','Region','RegionLabel'}));
    end
end

end