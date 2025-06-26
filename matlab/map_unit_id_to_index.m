
function index = map_unit_id_to_index(unit_id_to_index, unit_id)
    if isKey(unit_id_to_index, unit_id)
        index = unit_id_to_index(unit_id);
    else
        fields = split(unit_id, '|');
        if length(fields) == 5 && isKey(unit_id_to_index,[unit_id '||A'])
            index = unit_id_to_index([unit_id '||A']);
        elseif length(fields) == 8 && fields{8} == 'a'
            fields{8} = 'A';
            fix_id = strjoin(fields, '|');
            if isKey(unit_id_to_index, fix_id)
                index = unit_id_to_index(fix_id);
            end
        end
    end
end