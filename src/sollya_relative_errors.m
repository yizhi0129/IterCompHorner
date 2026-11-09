function rel_errors = sollya_relative_errors(p, x_value, values, precision)
%SOLLYA_RELATIVE_ERRORS Compute relative errors using Sollya arithmetic.

    if ~(isscalar(precision) && isfinite(precision) && precision >= 1 && ...
         precision == floor(precision))
        error('sollya_relative_errors:InvalidPrecision', ...
              'precision must be a positive integer.');
    end

    sollya_path = getenv('SOLLYA');
    if strlength(sollya_path) == 0
        sollya_path = '/opt/homebrew/bin/sollya';
    end

    script_file = [tempname, '.sollya'];
    cleanup = onCleanup(@() delete_if_exists(script_file));

    fid = fopen(script_file, 'w');
    if fid < 0
        error('sollya_relative_errors:FileOpenFailed', ...
              'Could not create temporary Sollya script.');
    end
    file_cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, 'prec = %d;\n', precision);
    fprintf(fid, 'display = decimal;\n');
    fprintf(fid, 'f = %s;\n', horner_expression(p));
    fprintf(fid, 'ref = evaluate(f, %.17g);\n', x_value);
    for k = 1:numel(values)
        fprintf(fid, 'v%d = %.17g;\n', k, values(k));
        fprintf(fid, ['rel%d = evaluate(abs(v%d - ref) / ', ...
                      'max(abs(ref), 2^(-1022)), 0);\n'], k, k);
        fprintf(fid, 'print("REL%d", rel%d);\n', k, k);
    end
    fprintf(fid, 'quit;\n');
    clear file_cleanup;

    command = sprintf('"%s" < "%s"', sollya_path, script_file);
    [status, output] = system(command);
    if status ~= 0
        error('sollya_relative_errors:SollyaFailed', ...
              'Sollya failed with status %d:\n%s', status, output);
    end

    rel_errors = NaN(size(values));
    for k = 1:numel(values)
        pattern = sprintf('REL%d\\s+([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?|[-+]?Inf|NaN)', k);
        token = regexp(output, pattern, 'tokens', 'once');
        if isempty(token)
            error('sollya_relative_errors:ParseFailed', ...
                  'Could not parse Sollya output for REL%d:\n%s', k, output);
        end
        rel_errors(k) = str2double(token{1});
    end
end

function expr = horner_expression(p)
    p = p(:).';
    expr = sprintf('%.17g', p(1));
    for k = 2:numel(p)
        expr = sprintf('((%s) * x + (%.17g))', expr, p(k));
    end
end

function delete_if_exists(filename)
    if exist(filename, 'file') == 2
        delete(filename);
    end
end
