% Generate richer analysis outputs for the custom full GECKO model.
%
% Outputs are written to:
%   GECKO-main/tutorials/full_ecModel/output/
%
% Generated files:
%   - full_model_summary.tsv
%   - full_model_flux_comparison.tsv
%   - full_model_ec_number_audit.tsv
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
rxnToEcPath = layout.rxnToEcCsv;
ecMapPath = fullfile(here, 'metaMet', 'data', 'raw', 'mapping_ec_numbers', 'mapping_ec_number_old_new.csv');
if ~isfile(fullPath) || ~isfile(kcatPath)
    error('Missing ecModel_full.yml or ecModel_full_kcat.yml. Run the full-model pipeline first.');
end

fullModel = readYAMLmodel(fullPath);
kcatModel = readYAMLmodel(kcatPath);
adapterParams = getFullAdapterParameters(layout);

fullSol = solveLP(fullModel);
kcatSol = solveLP(kcatModel);

expectedProtPoolLb = NaN;
if isfinite(adapterParams.sigma) && isfinite(adapterParams.Ptot) && isfinite(adapterParams.f)
    expectedProtPoolLb = -(adapterParams.sigma * adapterParams.Ptot * adapterParams.f * 1000);
end

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
    repmat(adapterParams.sigma, 2, 1), ...
    repmat(adapterParams.Ptot, 2, 1), ...
    repmat(adapterParams.f, 2, 1), ...
    repmat(adapterParams.gR_exp, 2, 1), ...
    repmat(adapterParams.bioRxn, 2, 1), ...
    repmat(expectedProtPoolLb, 2, 1), ...
    'VariableNames', {'model','n_rxns','n_mets','n_usage_prot','n_ec_rxns','n_nonzero_kcat','prot_pool_lb','fba_stat','objective','sigma','Ptot','f_param','gR_exp','bio_rxn','expected_prot_pool_lb'} ...
);
writetable(summaryTbl, fullfile(outDir, 'full_model_summary.tsv'), 'FileType', 'text', 'Delimiter', '\t');

commonRxns = intersect(string(fullModel.rxns), string(kcatModel.rxns), 'stable');
[~, fullIdx] = ismember(commonRxns, string(fullModel.rxns));
[~, kcatIdx] = ismember(commonRxns, string(kcatModel.rxns));
rxnToEcTbl = loadRxnToEcTable(rxnToEcPath);
ecMapInfo = loadEcOldNewInfo(ecMapPath);
[rxnBaseIds, ecNumbers, ecMapStatus, ecOldNewStatus] = annotateFluxRxns(commonRxns, rxnToEcTbl, ecMapInfo);
fluxTbl = table( ...
    commonRxns, ...
    rxnBaseIds, ...
    string(fullModel.rxnNames(fullIdx)), ...
    ecNumbers, ...
    ecMapStatus, ...
    ecOldNewStatus, ...
    fullSol.x(fullIdx), ...
    kcatSol.x(kcatIdx), ...
    abs(kcatSol.x(kcatIdx) - fullSol.x(fullIdx)), ...
    'VariableNames', {'rxn_id','rxn_base_id','rxn_name','ec_numbers','ec_map_status','ec_oldnew_status','flux_full','flux_full_kcat','abs_flux_delta'} ...
);
fluxTbl = sortrows(fluxTbl, 'abs_flux_delta', 'descend');
writetable(fluxTbl, fullfile(outDir, 'full_model_flux_comparison.tsv'), 'FileType', 'text', 'Delimiter', '\t');

ecAuditTbl = buildEcAuditTable(rxnToEcTbl, ecMapInfo);
writetable(ecAuditTbl, fullfile(outDir, 'full_model_ec_number_audit.tsv'), 'FileType', 'text', 'Delimiter', '\t');

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

disp(summaryTbl);
fprintf('GECKO full adapter params: sigma=%g, Ptot=%g, f=%g, gR_exp=%g, bioRxn=%s\n', ...
    adapterParams.sigma, adapterParams.Ptot, adapterParams.f, adapterParams.gR_exp, char(adapterParams.bioRxn));
fprintf('Expected prot_pool_exchange lower bound from GECKO params: %g\n', expectedProtPoolLb);
printStatusCounts(fluxTbl.ec_map_status, 'Flux-to-EC mapping status');
printStatusCounts(ecAuditTbl.ec_oldnew_status, 'rxn_to_ec old/new EC coverage');

fprintf('Wrote %s\n', fullfile(outDir, 'full_model_summary.tsv'));
fprintf('Wrote %s\n', fullfile(outDir, 'full_model_flux_comparison.tsv'));
fprintf('Wrote %s\n', fullfile(outDir, 'full_model_ec_number_audit.tsv'));
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

function params = getFullAdapterParameters(layout)
    params = struct('sigma', NaN, 'Ptot', NaN, 'f', NaN, 'gR_exp', NaN, 'bioRxn', "");
    try
        if exist('ModelAdapterManager', 'class') == 8
            adapterPath = fullfile(layout.matlabAdaptersDir, 'CarveMeFullModelAdapter.m');
            if isfile(adapterPath)
                ModelAdapterManager.setDefault(adapterPath, true);
                adapter = ModelAdapterManager.getDefault();
                p = adapter.getParameters();
                params.sigma = getNumericField(p, 'sigma');
                params.Ptot = getNumericField(p, 'Ptot');
                params.f = getNumericField(p, 'f');
                params.gR_exp = getNumericField(p, 'gR_exp');
                if isfield(p, 'bioRxn')
                    params.bioRxn = string(p.bioRxn);
                end
            end
        end
    catch ME
        warning('Could not read full-model adapter parameters for reporting: %s', ME.message);
    end
end

function value = getNumericField(s, fieldName)
    value = NaN;
    if isfield(s, fieldName)
        raw = s.(fieldName);
        if isnumeric(raw)
            value = double(raw);
        else
            value = str2double(string(raw));
        end
    end
end

function tbl = loadRxnToEcTable(path)
    if ~isfile(path)
        tbl = table(strings(0,1), strings(0,1), strings(0,1), 'VariableNames', {'rxn_id','ec_number','ec_raw'});
        return;
    end
    tbl = readtable(path, 'TextType', 'string');
    tbl.rxn_id = string(tbl.rxn_id);
    tbl.ec_number = strtrim(string(tbl.ec_number));
    if ismember('ec_raw', tbl.Properties.VariableNames)
        tbl.ec_raw = string(tbl.ec_raw);
    else
        tbl.ec_raw = repmat("", height(tbl), 1);
    end
end

function info = loadEcOldNewInfo(path)
    info = struct();
    info.oldVals = strings(0,1);
    info.newVals = strings(0,1);
    info.oldToNew = containers.Map('KeyType', 'char', 'ValueType', 'any');

    if ~isfile(path)
        return;
    end

    opts = detectImportOptions(path, 'FileType', 'text');
    opts = setvartype(opts, opts.VariableNames, 'string');
    opts.VariableNamingRule = 'preserve';
    tbl = readtable(path, opts);
    varNames = string(tbl.Properties.VariableNames);
    lowerNames = lower(varNames);

    oldIdx = find(contains(lowerNames, 'old'), 1, 'first');
    newIdx = find(contains(lowerNames, 'new'), 1, 'first');
    if isempty(oldIdx) || isempty(newIdx)
        return;
    end

    oldVals = strtrim(string(tbl.(varNames(oldIdx))));
    newVals = strtrim(string(tbl.(varNames(newIdx))));
    keep = oldVals ~= "" & newVals ~= "";
    oldVals = oldVals(keep);
    newVals = newVals(keep);
    info.oldVals = unique(oldVals);
    info.newVals = unique(newVals);

    uniqueOld = unique(oldVals);
    for i = 1:numel(uniqueOld)
        key = uniqueOld(i);
        replacements = unique(newVals(oldVals == key));
        info.oldToNew(char(key)) = replacements;
    end
end

function [baseIds, ecNumbers, mapStatus, oldNewStatus] = annotateFluxRxns(rxnIds, rxnToEcTbl, ecMapInfo)
    n = numel(rxnIds);
    baseIds = strings(n,1);
    ecNumbers = strings(n,1);
    mapStatus = strings(n,1);
    oldNewStatus = strings(n,1);

    for i = 1:n
        rxnId = string(rxnIds(i));
        baseId = normalizeFluxRxnId(rxnId);
        baseIds(i) = baseId;

        if startsWith(rxnId, 'usage_prot_')
            mapStatus(i) = "gecko_protein_usage";
            oldNewStatus(i) = "not_applicable";
            continue;
        elseif rxnId == "prot_pool_exchange"
            mapStatus(i) = "gecko_protein_pool";
            oldNewStatus(i) = "not_applicable";
            continue;
        elseif rxnId == "Growth"
            mapStatus(i) = "biomass_objective";
            oldNewStatus(i) = "not_applicable";
        end

        directMask = rxnToEcTbl.rxn_id == rxnId;
        baseMask = rxnToEcTbl.rxn_id == baseId;
        sourceMask = directMask;
        if any(directMask)
            mapStatus(i) = "mapped_direct";
        elseif any(baseMask)
            sourceMask = baseMask;
            mapStatus(i) = "mapped_by_base_id";
        elseif mapStatus(i) == ""
            mapStatus(i) = "no_ec_in_draft_model";
        end

        ecs = unique(strtrim(string(rxnToEcTbl.ec_number(sourceMask))));
        ecs(ecs == "") = [];
        if isempty(ecs)
            ecNumbers(i) = "";
            if oldNewStatus(i) == ""
                oldNewStatus(i) = "not_applicable";
            end
            continue;
        end

        ecNumbers(i) = join(ecs, ';');
        oldStatuses = classifyEcCoverage(ecs, ecMapInfo);
        if all(oldStatuses == "new_ec_in_map" | oldStatuses == "old_ec_in_map")
            oldNewStatus(i) = "all_listed_in_oldnew_map";
        elseif any(oldStatuses == "new_ec_in_map" | oldStatuses == "old_ec_in_map")
            oldNewStatus(i) = "partially_listed_in_oldnew_map";
        else
            oldNewStatus(i) = "not_listed_in_oldnew_map";
        end
    end
end

function tbl = buildEcAuditTable(rxnToEcTbl, ecMapInfo)
    if isempty(rxnToEcTbl)
        tbl = table(strings(0,1), strings(0,1), strings(0,1), strings(0,1), ...
            'VariableNames', {'rxn_id','ec_number','ec_oldnew_status','mapped_new_ec_numbers'});
        return;
    end

    tbl = unique(rxnToEcTbl(:, {'rxn_id','ec_number'}));
    n = height(tbl);
    status = strings(n,1);
    mappedNew = strings(n,1);
    for i = 1:n
        ec = strtrim(string(tbl.ec_number(i)));
        if ec == ""
            status(i) = "empty";
            mappedNew(i) = "";
            continue;
        end
        if any(ecMapInfo.oldVals == ec)
            status(i) = "old_ec_in_map";
            if isKey(ecMapInfo.oldToNew, char(ec))
                mappedNew(i) = join(string(ecMapInfo.oldToNew(char(ec))), ';');
            else
                mappedNew(i) = "";
            end
        elseif any(ecMapInfo.newVals == ec)
            status(i) = "new_ec_in_map";
            mappedNew(i) = ec;
        else
            status(i) = "not_listed_in_oldnew_map";
            mappedNew(i) = "";
        end
    end
    tbl.ec_oldnew_status = status;
    tbl.mapped_new_ec_numbers = mappedNew;
end

function baseId = normalizeFluxRxnId(rxnId)
    baseId = regexprep(string(rxnId), '_EXP_\d+$', '');
    baseId = regexprep(baseId, '_REV$', '');
end

function status = classifyEcCoverage(ecs, ecMapInfo)
    status = repmat("not_listed_in_oldnew_map", numel(ecs), 1);
    for i = 1:numel(ecs)
        ec = strtrim(string(ecs(i)));
        if any(ecMapInfo.oldVals == ec)
            status(i) = "old_ec_in_map";
        elseif any(ecMapInfo.newVals == ec)
            status(i) = "new_ec_in_map";
        end
    end
end

function printStatusCounts(values, label)
    values = string(values(:));
    values(values == "") = "(empty)";
    [u, ~, idx] = unique(values, 'stable');
    counts = accumarray(idx, 1);
    fprintf('%s:\n', label);
    for i = 1:numel(u)
        fprintf('  %s: %d\n', u(i), counts(i));
    end
end
