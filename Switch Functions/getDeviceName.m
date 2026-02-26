function name = getDeviceName(T, ii, tech, fallbackIdx)
    if istable(T) && ~isempty(T) && ismember('Name', T.Properties.VariableNames)
        try
            nm = string(T.Name(ii));
            if strlength(nm) == 0 || nm == "missing"
                nm = "device#" + string(fallbackIdx);
            end
        catch
            nm = "device#" + string(fallbackIdx);
        end
    else
        nm = "device#" + string(fallbackIdx);
    end
    name = tech + ": " + nm;
end