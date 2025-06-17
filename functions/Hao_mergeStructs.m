function mergedStruct = Hao_mergeStructs(baseStruct, newStruct)
    % Ensure inputs are both structures
    if ~isstruct(baseStruct) || ~isstruct(newStruct)
        error('Both inputs must be structures.');
    end
    
    % Get field names from newStruct
    fields = fieldnames(newStruct);
    
    % Merge fields into baseStruct
    for i = 1:numel(fields)
        baseStruct.(fields{i}) = newStruct.(fields{i});
    end
    
    mergedStruct = baseStruct;
end
