% Generate richer analysis outputs for the custom full GECKO model.
%
% Outputs are written to:
%   GECKO-main/tutorials/full_ecModel/output/
%
% Generated files:
%   - full_model_summary.tsv
%   - full_model_flux_comparison.tsv
%   - full_model_kcat_distribution.tsv
%   - full_model_top_enzyme_usage.tsv
%   - full_model_overview.pdf

clearvars;

here = resolve_workspace_root();
layout = resolve_modelyml_layout(here);
geckoRoot = resolve_gecko_root(here);
if isfolder(geckoRoot)
    addpath(genpath(fullfile(geckoRoot, 'src')));
end

if exist('readYAMLmodel', 'file') ~= 2
    error('readYAMLmodel not found on MATLAB path.');
end
if exist('solveLP', 'file') ~= 2
    error('solveLP not found on MATLAB path.');
end

outDir = fullfile(geckoRoot, 'tutorials', 'full_ecModel', 'output');
if ~isfolder(outDir)
    mkdir(outDir);
end

fullPath = layout.fullEcModel;
kcatPath = layout.fullEcModelKcat;
if ~isfile(fullPath) || ~isfile(kcatPath)
    error('Missing ecModel_full.yml or ecModel_full_kcat.yml. Run the full-model pipeline first.');
end

fullModel = readYAMLmodel(fullPath);
kcatModel = readYAMLmodel(kcatPath);

fullSol = solveLP(fullModel);
kcatSol = solveLP(kcatModel);

summaryTbl = table( ...
    string({'ecModel_full'; 'ecModel_full_kcat'}), ...
    [numel(fullModel.rxns); numel(kcatModel.rxns)], ...
    [numel(fullModel.mets); numel(kcatModel.mets)], ...
    [sum(startsWith(string(fullModel.rxns), 'usage_prot_')); sum(startsWith(string(kcatModel.rxns), 'usage_prot_'))], ...
    [countEcRxns(fullModel); countEcRxns(kcatModel)], ...
    [countNonzeroKcats(fullModel); countNonzeroKcats(kcatModel)], ...
    [getProtPoolLb(fullModel); getProtPoolLb(kcatModel)], ...
    [fullSol.stat; kcatSol.stat], ...
    [fullSol.f; kcatSol.f], ...
    'VariableNames', {'model','n_rxns','n_mets','n_usage_prot','n_ec_rxns','n_nonzero_kcat','prot_pool_lb','fba_stat','objective'} ...
);
writetable(summaryTbl, fullfile(outDir, 'full_model_summary.tsv'), 'FileType', 'text', 'Delimiter', '\t');

commonRxns = intersect(string(fullModel.rxns), string(kcatModel.rxns), 'stable');
[~, fullIdx] = ismember(commonRxns, string(fullModel.rxns));
[~, kcatIdx] = ismember(commonRxns, string(kcatModel.rxns));
fluxTbl = table( ...
    commonRxns, ...
    string(fullModel.rxnNames(fullIdx)), ...
    fullSol.x(fullIdx), ...
    kcatSol.x(kcatIdx), ...
    abs(kcatSol.x(kcatIdx) - fullSol.x(fullIdx)), ...
    'VariableNames', {'rxn_id','rxn_name','flux_full','flux_full_kcat','abs_flux_delta'} ...
);
fluxTbl = sortrows(fluxTbl, 'abs_flux_delta', 'descend');
writetable(fluxTbl, fullfile(outDir, 'full_model_flux_comparison.tsv'), 'FileType', 'text', 'Delimiter', '\t');

if isfield(kcatModel, 'ec') && isfield(kcatModel.ec, 'kcat')
    kcatVals = kcatModel.ec.kcat(:);
    kcatVals = kcatVals(kcatVals > 0 & isfinite(kcatVals));
else
    kcatVals = zeros(0,1);
end
if isempty(kcatVals)
    kcatTbl = table([], [], 'VariableNames', {'kcat', 'log10_kcat'});
else
    kcatTbl = table(kcatVals, log10(kcatVals), 'VariableNames', {'kcat', 'log10_kcat'});
end
writetable(kcatTbl, fullfile(outDir, 'full_model_kcat_distribution.tsv'), 'FileType', 'text', 'Delimiter', '\t');

enzymeTbl = table();
if exist('enzymeUsage', 'file') == 2 && exist('reportEnzymeUsage', 'file') == 2
    try
        usageData = enzymeUsage(kcatModel, kcatSol.x);
        usageReport = reportEnzymeUsage(kcatModel, usageData, 0.9, min(15, numel(usageData.absUsage)));
        enzymeTbl = usageReport.topAbsUsage;
    catch ME
        warning('Skipping enzyme usage report: %s', ME.message);
    end
end
if isempty(enzymeTbl)
    enzymeTbl = table(strings(0,1), strings(0,1), zeros(0,1), zeros(0,1), ...
        'VariableNames', {'protID','geneID','absUsage','percUsage'});
end
writetable(enzymeTbl, fullfile(outDir, 'full_model_top_enzyme_usage.tsv'), 'FileType', 'text', 'Delimiter', '\t');

fig = figure('Visible', 'off', 'Position', [100 100 1400 900]);
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
bar(categorical(summaryTbl.model), summaryTbl.objective);
ylabel('Objective value');
title('Growth / objective comparison');

nexttile;
plotFluxScatter(fullSol.x(fullIdx), kcatSol.x(kcatIdx));
title('Full vs kcat full fluxes');

nexttile;
if isempty(kcatVals)
    text(0.5, 0.5, 'No nonzero kcat values', 'HorizontalAlignment', 'center');
    axis off;
else
    histogram(log10(kcatVals), 30);
    xlabel('log_{10}(kcat)');
    ylabel('Count');
    title('Nonzero kcat distribution');
end

nexttile;
plotTopEnzymeUsage(enzymeTbl);
title('Top enzyme usage in kcat model');

exportgraphics(fig, fullfile(outDir, 'full_model_overview.pdf'), 'ContentType', 'vector');
close(fig);

fprintf('Wrote %s\n', fullfile(outDir, 'full_model_summary.tsv'));
fprintf('Wrote %s\n', fullfile(outDir, 'full_model_flux_comparison.tsv'));
fprintf('Wrote %s\n', fullfile(outDir, 'full_model_kcat_distribution.tsv'));
fprintf('Wrote %s\n', fullfile(outDir, 'full_model_top_enzyme_usage.tsv'));
fprintf('Wrote %s\n', fullfile(outDir, 'full_model_overview.pdf'));

function n = countEcRxns(model)
    n = 0;
    if isfield(model, 'ec') && isfield(model.ec, 'rxns')
        n = numel(model.ec.rxns);
    end
end

function n = countNonzeroKcats(model)
    n = 0;
    if isfield(model, 'ec') && isfield(model.ec, 'kcat')
        k = model.ec.kcat;
        n = sum(k > 0);
    end
end

function v = getProtPoolLb(model)
    v = NaN;
    idx = strcmp(string(model.rxns), 'prot_pool_exchange');
    if any(idx)
        v = model.lb(find(idx,1));
    end
end

function plotFluxScatter(x, y)
    x = abs(x);
    y = abs(y);
    x(x < 1e-12) = 1e-12;
    y(y < 1e-12) = 1e-12;
    loglog(x, y, '.', 'MarkerSize', 8);
    hold on;
    lims = [min([x; y]), max([x; y])];
    plot(lims, lims, 'k--');
    hold off;
    xlabel('|flux| full');
    ylabel('|flux| full+kcat');
    grid on;
end

function plotTopEnzymeUsage(enzymeTbl)
    if isempty(enzymeTbl) || height(enzymeTbl) == 0 || ~ismember('protID', enzymeTbl.Properties.VariableNames)
        text(0.5, 0.5, 'No enzyme usage report available', 'HorizontalAlignment', 'center');
        axis off;
        return;
    end
    keep = isfinite(enzymeTbl.absUsage);
    enzymeTbl = enzymeTbl(keep,:);
    if height(enzymeTbl) == 0
        text(0.5, 0.5, 'No enzyme usage report available', 'HorizontalAlignment', 'center');
        axis off;
        return;
    end
    n = min(10, height(enzymeTbl));
    enzymeTbl = enzymeTbl(1:n,:);
    barh(enzymeTbl.absUsage(end:-1:1));
    set(gca, 'YTick', 1:n, 'YTickLabel', cellstr(enzymeTbl.protID(end:-1:1)));
    xlabel('Absolute usage');
    ylabel('Protein');
end
