function e_vectors = accsumk_vector_full(vectors)
    nvec = numel(vectors);
    if nvec == 0
        e_vectors = {0};
        return;
    end

    lengths = zeros(1, nvec);
    for k = 1:nvec
        v = vectors{k};
        if isempty(v)
            vectors{k} = 0;
            lengths(k) = 1;
        else
            vectors{k} = v(:).';
            lengths(k) = numel(vectors{k});
        end
    end

    max_len = max(lengths);
    aligned = zeros(nvec, max_len);
    for k = 1:nvec
        v = vectors{k};
        aligned(k, end - numel(v) + 1:end) = v;
    end

    expanded = zeros(nvec, max_len);
    for j = 1:max_len
        col = aligned(:, j).';
        if any(col ~= 0)
            expanded(:, j) = AccSumK(col, nvec).';
        end
    end

    e_vectors = cell(1, nvec);
    for k = 1:nvec
        e_vectors{k} = expanded(k, :);
    end
end
