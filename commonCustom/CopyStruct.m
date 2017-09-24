function Scopy = CopyStruct(S)
    % Drew Chap 6/5/17
    % forces a copy of a structure (otherwise MATLAB does lazy copy)
    SVfields = fieldnames(S);
    for f = 1:numel(SVfields)
        Scopy.(SVfields{f}) = 1*S.(SVfields{f});  % multiply by 1 to force a copy
    end
end