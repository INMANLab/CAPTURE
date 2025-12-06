function [descriptionAll, descriptionUnique] = ExtractDescriptions(descriptionCol,pattern)
    descriptionAll = regexp(descriptionCol, pattern, 'tokens');
    
    descriptionUnique = descriptionAll;
    while(iscell(descriptionUnique))
        descriptionUnique = horzcat(descriptionUnique{:});
    end
    descriptionUnique = unique(descriptionUnique)';

    descriptionAll = cellfun(@(x) combineStrings(x),descriptionAll,'UniformOutput',false);
    descriptionAll = vertcat(descriptionAll{:});
end

function s = combineStrings(s)
    if(iscell(s))
        s = cellfun(@(y) ("#"+y+";"),s,'UniformOutput',false);
        s = strjoin(horzcat(s{:}));
    elseif(isempty(s))
        s = "";
    else
        s = string(s);
    end
end