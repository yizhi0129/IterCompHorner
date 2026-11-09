function plot_iter_CompHorner_sweeps()
%PLOT_ITER_COMPHORNER_SWEEPS Export four plots from the saved sweep reports.

    root = fileparts(fileparts(mfilename('fullpath')));
    results_dir = fullfile(root, 'results');
    fig_dir = fullfile(root, 'fig');
    if ~isfolder(fig_dir)
        mkdir(fig_dir);
    end

    prefix = 'iter_CompHorner_gamma0_mex_degree';
    suffix = '_sweep_1e2_1e35_relerr_output.txt';
    random_data = read_sweep(fullfile(results_dir, [prefix '50_random_x_v' suffix]), 2);
    degrees = [5, 10, 50, 100, 300, 500];
    fixed_data = cell(size(degrees));
    for k = 1:numel(degrees)
        filename = [prefix num2str(degrees(k)) suffix];
        fixed_data{k} = sortrows(read_sweep(fullfile(results_dir, filename), 0), 1);
    end

    save_plot({random_data}, [], 2, fullfile(fig_dir, 'random_x_v_iterations.png'));
    save_plot({random_data}, [], 3, fullfile(fig_dir, 'random_x_v_relative_error.png'));
    save_plot(fixed_data, degrees, 2, fullfile(fig_dir, 'degree_sweep_iterations.png'));
    save_plot(fixed_data, degrees, 3, fullfile(fig_dir, 'degree_sweep_relative_error.png'));
end

function data = read_sweep(filename, offset)
    lines = strip(readlines(filename));
    pattern = ['^' repmat('\S+\s+', 1, offset) '1e\d+\s+'];
    matches = regexp(cellstr(lines), pattern, 'once');
    rows = lines(~cellfun(@isempty, matches));
    assert(~isempty(rows), 'No sweep data found in %s.', filename);
    columns = split(regexprep(rows, '\s+', ' '));
    data = str2double(columns(:, offset + [2, 4, 6]));
end

function save_plot(series, degrees, column, filename)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 520]);
    cleanup = onCleanup(@() close(fig)); %#ok<NASGU>
    ax = axes(fig);
    ax.FontSize = 32;
    ax.LabelFontSizeMultiplier = 1.1;
    ax.TitleFontSizeMultiplier = 1.2;
    hold(ax, 'on');
    is_error = column == 3;
    if isempty(degrees)
        title(ax, 'Degree 50, random x and v');
    else
        title(ax, 'Fixed x and v, by polynomial degree');
    end

    for k = 1:numel(series)
        x = series{k}(:, 1);
        y = series{k}(:, column);
        y(~isfinite(x) | x <= 0) = NaN;
        if is_error
            y(y <= 0) = NaN;
        end
        if isempty(degrees)
            scatter(ax, x, y, 62, [0.05 0.25 0.51], 'filled', 'MarkerFaceAlpha', 0.82);
        else
            plot(ax, x, y, '-o', 'LineWidth', 3, 'MarkerSize', 3.5, ...
                'DisplayName', sprintf('Degree %d', degrees(k)));
        end
    end

    ax.XScale = 'log';
    xlabel(ax, 'Actual condition number');
    values = vertcat(series{:});
    valid = isfinite(values(:, 1)) & values(:, 1) > 0 & isfinite(values(:, column));
    plotted_y = values(valid, column);
    if is_error
        ax.YScale = 'log';
        ylabel(ax, 'Relative error');
        u = 2^-53;
        positive_y = [plotted_y(plotted_y > 0); u];
        ylim(ax, [0.8 * min(positive_y), 1.6 * max(positive_y)]);
        yline(ax, u, '--k', 'u = 2^{-53}', 'LineWidth', 2, ...
            'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
    else
        ylabel(ax, 'Iteration count');
        if ~isempty(plotted_y)
            yticks(ax, min(plotted_y):max(plotted_y));
            ylim(ax, [min(plotted_y) - 0.3, max(plotted_y) + 0.3]);
        end
    end
    if ~isempty(degrees)
        legend(ax, 'Location', 'eastoutside', 'FontSize', 30);
    end
    grid(ax, 'on');
    box(ax, 'on');
    exportgraphics(ax, filename, 'Resolution', 300);
end
