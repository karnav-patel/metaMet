% Apply metaMet kcats to the full GECKO ecModel and finalize the protein-
% constrained full model so it is ready for simulation.
%
% Inputs:
%   - ecModel_full.yml
%   - rxn_to_ec.csv
%   - metaMet/data/processed/overview/kcat_aggregate.csv
%   - gecko_full_adapter/data/uniprot.tsv
%   - gecko_full_adapter/data/pseudoRxns.tsv
%
% Outputs:
%   - gecko_full_adapter/data/customKcats.tsv
%   - ecModel_full_kcat.yml

clearvars;

workspaceRoot = resolve_workspace_root();
geckoRoot = resolve_gecko_root(workspaceRoot);
layout = resolve_modelyml_layout(workspaceRoot);
addpath(genpath(fullfile(geckoRoot, 'src')));
addpath(workspaceRoot);

ModelAdapterManager.setDefault(fullfile(layout.matlabAdaptersDir, 'CarveMeFullModelAdapter.m'), true);
run('prepare_gecko_full_adapter.m');
adapter = ModelAdapterManager.getDefault();
params = adapter.getParameters();

inModelPath = layout.fullEcModel;
if exist('readYAMLmodel', 'file') ~= 2 || ~isfile(inModelPath)
    error('Missing ecModel_full.yml or readYAMLmodel not on path.');
end
model = readYAMLmodel(inModelPath);

if exist('getECfromGEM', 'file') == 2
    model = getECfromGEM(model);
end
if exist('getECfromDatabase', 'file') == 2
    try
        noEc = ~isfield(model.ec, 'eccodes') || isempty(model.ec.eccodes);
        if islogical(noEc) && noEc
            model = getECfromDatabase(model, [], 'ignore', adapter);
        else
            missingEc = cellfun(@isempty, model.ec.eccodes);
            if any(missingEc)
                model = getECfromDatabase(model, missingEc, 'ignore', adapter);
            end
        end
    catch ME
        warning('Database-based EC assignment skipped: %s', ME.message);
    end
end

rxnToEcPath = layout.rxnToEcCsv;
kcatAggPath = fullfile(workspaceRoot, 'metaMet', 'data', 'processed', 'overview', 'kcat_aggregate.csv');
if ~isfile(rxnToEcPath)
    error('Missing %s', rxnToEcPath);
end
if ~isfile(kcatAggPath)
    error('Missing %s', kcatAggPath);
end

Tmap = readtable(rxnToEcPath, 'TextType', 'string');
Tagg = readtable(kcatAggPath, 'TextType', 'string');

ecMapPath = fullfile(workspaceRoot, 'metaMet', 'data', 'raw', 'mapping_ec_numbers', 'mapping_ec_number_old_new.csv');
if isfile(ecMapPath)
    ecMap = loadEcOldNewMap(ecMapPath);
    fprintf('Loaded EC old->new map (%d old ECs).\n', ecMap.Count);
else
    ecMap = [];
end

[ecKeys, ecKcats] = explodeEcToKcat(Tagg.ec_number, Tagg.kcat_final, ecMap);
if isempty(ecKeys)
    error('No valid EC->kcat pairs found in %s (after normalization).', kcatAggPath);
end

[uEc, ~, idx] = unique(ecKeys);
maxK = accumarray(idx, ecKcats, [], @max);
ecToK = containers.Map(cellstr(uEc), num2cell(maxK));

rxns = string(Tmap.rxn_id);
ecsRaw = string(Tmap.ec_number);
rxnList = unique(rxns);
rxnK = zeros(numel(rxnList), 1);

for i = 1:numel(rxnList)
    r = rxnList(i);
    theseRaw = ecsRaw(rxns == r);
    these = normalizeEcList(theseRaw, ecMap);
    kVals = [];
    for j = 1:numel(these)
        e = char(these(j));
        if isKey(ecToK, e)
            kVals(end+1,1) = ecToK(e); %#ok<AGROW>
        end
    end
    if ~isempty(kVals)
        rxnK(i) = max(kVals);
    end
end

hasK = rxnK > 0;
fprintf('Found kcats for %d/%d reactions with ECs.\n', sum(hasK), numel(rxnList));

customPath = fullfile(params.path, 'data', 'customKcats.tsv');
if ~isfolder(fullfile(params.path, 'data'))
    mkdir(fullfile(params.path, 'data'));
end

ecRxnNoExp = regexprep(string(model.ec.rxns), '_EXP_\d+$', '');
ecRxnBases = unique(ecRxnNoExp);

fid = fopen(customPath, 'wt');
if fid < 0
    error('Could not write %s', customPath);
end
fprintf(fid, 'proteins\tgenes\tgene_name\tkcat\trxns\tnotes\tstoicho\n');
for i = 1:numel(rxnList)
    if ~hasK(i)
        continue;
    end
    r = char(rxnList(i));
    k = rxnK(i);
    rxnField = r;
    if any(ecRxnBases == (string(r) + "_REV"))
        rxnField = [r ',' r '_REV'];
    end
    fprintf(fid, '\t\t\t%.6g\t%s\tmetamet_kcat_aggregate\t\n', k, rxnField);
end
fclose(fid);
fprintf('Wrote %s\n', customPath);

if exist('applyCustomKcats', 'file') ~= 2
    error('applyCustomKcats not found on path.');
end
[model, rxnUpdated] = applyCustomKcats(model, customPath, adapter);
fprintf('Updated %d ec reactions with custom kcats.\n', numel(rxnUpdated));

if exist('getKcatAcrossIsozymes', 'file') == 2
    model = getKcatAcrossIsozymes(model);
end
if exist('getStandardKcat', 'file') == 2
    [model, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(model, adapter);
    fprintf('Assigned standard kcat %.6g (MW %.6g) to %d missing-GPR reactions.\n', standardKcat, standardMW, numel(rxnsMissingGPR));
end
if exist('applyKcatConstraints', 'file') == 2
    model = applyKcatConstraints(model);
end
if exist('setProtPoolSize', 'file') == 2
    model = setProtPoolSize(model, [], [], [], adapter);
end

model = ensureNestedSubsystems(model);
outPath = layout.fullEcModelKcat;
writeYAMLmodel(model, outPath);
fprintf('Wrote %s\n', outPath);

metaMetModelsDir = layout.metaMetModelsDir;
if ~isfolder(metaMetModelsDir)
    mkdir(metaMetModelsDir);
end
copyfile(outPath, fullfile(metaMetModelsDir, 'ecModel_full_kcat.yml'));
copyfile(inModelPath, fullfile(metaMetModelsDir, 'ecModel_full.yml'));
copyfile(customPath, fullfile(metaMetModelsDir, 'customKcats_full.tsv'));
fprintf('Copied full-model outputs into %s\n', metaMetModelsDir);

tutorialOutputDir = fullfile(geckoRoot, 'tutorials', 'full_ecModel', 'output');
if isfolder(tutorialOutputDir)
    copyfile(inModelPath, fullfile(tutorialOutputDir, 'ecModel_full.yml'));
    copyfile(outPath, fullfile(tutorialOutputDir, 'ecModel_full_kcat.yml'));
    copyfile(customPath, fullfile(tutorialOutputDir, 'customKcats_full.tsv'));
    if isfile(rxnToEcPath)
        copyfile(rxnToEcPath, fullfile(tutorialOutputDir, 'rxn_to_ec_full.csv'));
    end
    fprintf('Copied full-model outputs into %s\n', tutorialOutputDir);
end

function [ecKeys, ecKcats] = explodeEcToKcat(ecCol, kcatCol, ecMap)
    if nargin < 3
        ecMap = [];
    end
    ecCol = string(ecCol);
    kcatCol = string(kcatCol);
    kcatNum = str2double(kcatCol);
    ecKeys = strings(0, 1);
    ecKcats = zeros(0, 1);
    for i = 1:numel(ecCol)
        if isnan(kcatNum(i))
            continue;
        end
        s = normalizeText(ecCol(i));
        if s == ""
            continue;
        end
        tokens = extractEcTokens(s);
        ecs = normalizeEcTokens(tokens, ecMap);
        if isempty(ecs)
            continue;
        end
        ecKeys = [ecKeys; ecs(:)]; %#ok<AGROW>
        ecKcats = [ecKcats; repmat(kcatNum(i), numel(ecs), 1)]; %#ok<AGROW>
    end
    keep = ecKeys ~= "" & ~ismissing(ecKeys) & ~isnan(ecKcats);
    ecKeys = ecKeys(keep);
    ecKcats = ecKcats(keep);
end

function ecs = normalizeEcList(raw, ecMap)
    if nargin < 2
        ecMap = [];
    end
    raw = string(raw);
    ecs = strings(0, 1);
    for i = 1:numel(raw)
        s = normalizeText(raw(i));
        if s == ""
            continue;
        end
        tokens = extractEcTokens(s);
        found = normalizeEcTokens(tokens, ecMap);
        if ~isempty(found)
            ecs = [ecs; found(:)]; %#ok<AGROW>
        end
    end
    ecs = unique(ecs);
end

function s = normalizeText(s)
    s = string(s);
    if ismissing(s)
        s = "";
        return;
    end
    s = strtrim(s);
    if s == ""
        return;
    end
    sLower = lower(s);
    if sLower == "nan" || sLower == "na" || sLower == "none" || sLower == "null" || sLower == "<missing>"
        s = "";
        return;
    end
    s = erase(s, '"');
    s = erase(s, "'");
    s = strtrim(s);
end

function tokens = extractEcTokens(s)
    s = char(string(s));
    m = regexp(s, '(?<!\d)(\d+)\.(\d+)\.(\d+)\.([A-Za-z0-9]+)(?!\d)', 'match');
    if isempty(m)
        tokens = strings(0, 1);
    else
        tokens = unique(string(m(:)));
    end
end

function ecs = normalizeEcTokens(tokens, ecMap)
    tokens = string(tokens(:));
    out = strings(0, 1);
    for i = 1:numel(tokens)
        t = normalizeText(tokens(i));
        if t == ""
            continue;
        end
        if ~isempty(ecMap) && isKey(ecMap, char(t))
            t = string(ecMap(char(t)));
        end
        if ~isempty(regexp(char(t), '^\d+\.\d+\.\d+\.\d+$', 'once'))
            out(end+1,1) = t; %#ok<AGROW>
        end
    end
    ecs = unique(out);
end

function ecMap = loadEcOldNewMap(ecMapPath)
    T = readtable(ecMapPath, 'TextType', 'string');
    oldCol = string(T{:,1});
    newCol = string(T{:,2});
    keep = oldCol ~= "" & newCol ~= "" & ~ismissing(oldCol) & ~ismissing(newCol);
    ecMap = containers.Map(cellstr(oldCol(keep)), cellstr(newCol(keep)));
end

function model = ensureNestedSubsystems(model)
    if ~isfield(model, 'rxns')
        return;
    end
    nRxns = numel(model.rxns);
    if ~isfield(model, 'subSystems') || isempty(model.subSystems)
        model.subSystems = repmat({{''}}, nRxns, 1);
        return;
    end
    ss = model.subSystems;
    if isstring(ss)
        ss = cellstr(ss);
    elseif ischar(ss)
        ss = {ss};
    end
    if ~iscell(ss)
        ss = cellstr(string(ss));
    end
    ss = ss(:);
    out = repmat({{''}}, nRxns, 1);
    for i = 1:nRxns
        if i > numel(ss)
            out{i} = {''};
            continue;
        end
        entry = ss{i};
        if isstring(entry)
            entry = cellstr(entry);
        elseif ischar(entry)
            entry = {entry};
        elseif ~iscell(entry)
            entry = {char(string(entry))};
        end
        entry = entry(:)';
        entry = entry(~cellfun(@(x) isempty(strtrim(char(string(x)))), entry));
        if isempty(entry)
            entry = {''};
        else
            entry = cellfun(@(x) char(string(x)), entry, 'UniformOutput', false);
        end
        out{i} = entry;
    end
    model.subSystems = out;
end