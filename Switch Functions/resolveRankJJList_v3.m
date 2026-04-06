%% Match with v3 for analyzePriSwitches_v3.
function jj_list_out = resolveRankJJList_v3(mode, jj_list_in, max_para, selected_para)
% Supports:
%   1) scalar-jj mode:
%        jj_list_out = [1 2 3 ...]
%
%   2) pair-jj mode:
%        jj_list_out = [jj25 jj75;
%                       jj25 jj75;
%                       ...]
%
% Rules:
%   mode == 1:
%       - scalar mode: keep entries within [1, max_para]
%       - pair mode:   keep rows whose two entries are within [1, max_para]
%                      if jj_list_in is empty, generate full grid
%
%   mode == 2:
%       - scalar mode: restrict to selected_para
%       - pair mode:   restrict to rows in selected_para
%
% Notes:
%   - selected_para can be either:
%         vector      -> scalar mode
%         N x 2 matrix -> pair mode
%
%   - jj_list_in can be empty; then the function falls back sensibly.

    %% ---------------- Detect scalar vs pair mode ----------------
    isPairMode = ~isempty(selected_para) && ismatrix(selected_para) && size(selected_para,2) == 2;

    % If selected_para is empty, infer from jj_list_in if possible
    if isempty(selected_para)
        if ~isempty(jj_list_in) && ismatrix(jj_list_in) && size(jj_list_in,2) == 2
            isPairMode = true;
        else
            isPairMode = false;
        end
    end

    %% ---------------- Normalize input ----------------
    if isPairMode
        if isempty(jj_list_in)
            jj_in = zeros(0,2);
        else
            if size(jj_list_in,2) ~= 2
                error('In pair mode, jj_list_in must be an N x 2 matrix.');
            end
            jj_in = unique(jj_list_in, 'rows', 'stable');
        end
    else
        jj_in = unique(jj_list_in(:).', 'stable');
    end

    %% ---------------- Resolve by mode ----------------
    if mode == 1

        if isPairMode
            if isempty(jj_in)
                % default full grid
                [J25, J75] = ndgrid(1:max_para, 1:max_para);
                jj_list_out = [J25(:), J75(:)];
            else
                valid = jj_in(:,1) >= 1 & jj_in(:,1) <= max_para & ...
                        jj_in(:,2) >= 1 & jj_in(:,2) <= max_para;
                jj_list_out = jj_in(valid, :);

                if isempty(jj_list_out)
                    [J25, J75] = ndgrid(1:max_para, 1:max_para);
                    jj_list_out = [J25(:), J75(:)];
                end
            end

        else
            jj_list_out = jj_in(jj_in >= 1 & jj_in <= max_para);
            if isempty(jj_list_out)
                jj_list_out = 1:max_para;
            end
        end

    elseif mode == 2

        if isPairMode
            if isempty(selected_para) || size(selected_para,2) ~= 2
                error('In pair mode and mode==2, selected_para must be an N x 2 matrix.');
            end

            sel = unique(selected_para, 'rows', 'stable');

            if isempty(jj_in)
                jj_list_out = sel;
            else
                tf = ismember(jj_in, sel, 'rows');
                overlap = jj_in(tf, :);

                if isempty(overlap)
                    jj_list_out = sel;
                else
                    jj_list_out = overlap;
                end
            end

        else
            sel = unique(selected_para(:).', 'stable');
            overlap = jj_in(ismember(jj_in, sel));

            if isempty(overlap)
                jj_list_out = sel;
            else
                jj_list_out = overlap;
            end
        end

    else
        error('mode must be 1 or 2');
    end
end