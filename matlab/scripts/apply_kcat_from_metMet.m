% Apply kcats to the GECKO ecModel using metaMet aggregated kcat table.
%
% Inputs:
%   - ecModel.yml (from build_ecmodel_gecko.m)
%   - rxn_to_ec.csv (exported earlier)
%   - metaMet/data/processed/overview/kcat_aggregate.csv
%
% Outputs:
%   - gecko_adapter/data/customKcats.tsv
%   - ecModel_kcat.yml (updated ecModel with ec.kcat filled + constraints applied)

clearvars;

workspaceRoot = resolve_workspace_root();
geckoRoot = resolve_gecko_root(workspaceRoot);
layout = resolve_modelyml_layout(workspaceRoot);

% Add GECKO to path
addpath(genpath(fullfile(geckoRoot, 'src')));

% Ensure adapter + UniProt-like DB exist
ModelAdapterManager.setDefault(fullfile(layout.matlabAdaptersDir, 'CarveMeModelAdapter.m'), true);
run('prepare_gecko_adapter.m');
adapter = ModelAdapterManager.getDefault();
params = adapter.getParameters();

% Load ecModel
inModelPath = layout.lightEcModel;
if exist('readYAMLmodel', 'file') ~= 2 || ~isfile(inModelPath)
    error('Missing ecModel.yml or readYAMLmodel not on path.');
end
model = readYAMLmodel(inModelPath);

% Populate ecModel.ec.eccodes from ecModel.eccodes (copied from GEM)
if exist('getECfromGEM', 'file') == 2
    model = getECfromGEM(model);
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

% Optional: normalize retired/old EC numbers using metaMet mapping
ecMapPath = fullfile(workspaceRoot, 'metaMet', 'data', 'raw', 'mapping_ec_numbers', 'mapping_ec_number_old_new.csv');
if isfile(ecMapPath)
    ecMap = loadEcOldNewMap(ecMapPath);
    fprintf('Loaded EC old->new map (%d old ECs).\n', ecMap.Count);
else
    ecMap = [];
end

% Build EC -> kcat map from Tagg.kcat_final.
% Robustly handle missing values, quoted strings, and multi-EC rows by extracting
% valid EC patterns and expanding them into one row per EC.
[ecKeys, ecKcats] = explodeEcToKcat(Tagg.ec_number, Tagg.kcat_final, ecMap);
if isempty(ecKeys)
    error('No valid EC->kcat pairs found in %s (after normalization).', kcatAggPath);
end

% If duplicate ECs exist, keep the maximum (conservative for growth feasibility)
[uEc, ~, idx] = unique(ecKeys);
maxK = accumarray(idx, ecKcats, [], @max);
ecToK = containers.Map(cellstr(uEc), num2cell(maxK));

% Compute reaction -> kcat as max over its ECs
rxns = string(Tmap.rxn_id);
ecsRaw  = string(Tmap.ec_number);

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
    else
        rxnK(i) = 0;
    end
end

hasK = rxnK > 0;
fprintf('Found kcats for %d/%d reactions with ECs.\n', sum(hasK), numel(rxnList));

% Build customKcats.tsv (rxn-only entries) for GECKO applyCustomKcats
customPath = fullfile(params.path, 'data', 'customKcats.tsv');
if ~isfolder(fullfile(params.path, 'data'))
    mkdir(fullfile(params.path, 'data'));
end

% Determine which reactions have a reversible version in ecModel (so we can also set *_REV)
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

% Apply to model
if exist('applyCustomKcats', 'file') ~= 2
    error('applyCustomKcats not found on path.');
end
[model, rxnUpdated] = applyCustomKcats(model, customPath, adapter);
fprintf('Updated %d ec reactions with custom kcats.\n', numel(rxnUpdated));

% Export updated ecModel
model = ensureNestedSubsystems(model);
outPath = layout.lightEcModelKcat;
writeYAMLmodel(model, outPath);
fprintf('Wrote %s\n', outPath);

% Also copy outputs into the metaMet project tree (so everything is co-located)
metaMetModelsDir = layout.metaMetModelsDir;
if ~isfolder(metaMetModelsDir)
    mkdir(metaMetModelsDir);
end
copyfile(outPath, fullfile(metaMetModelsDir, 'ecModel_kcat.yml'));
copyfile(inModelPath, fullfile(metaMetModelsDir, 'ecModel.yml'));
copyfile(layout.draftModelYaml, fullfile(metaMetModelsDir, 'model.yml'));
copyfile(customPath, fullfile(metaMetModelsDir, 'customKcats.tsv'));
fprintf('Copied models into %s\n', metaMetModelsDir);

function [ecKeys, ecKcats] = explodeEcToKcat(ecCol, kcatCol, ecMap)
% Extract valid EC patterns (a.b.c.d, digits only) from arbitrary strings.
% Expands multi-EC rows into one EC per key.
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
% Given a string array (potentially with multi-EC entries), return a unique
% list of extracted EC patterns.
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
% Normalize arbitrary text field to a safe string for further parsing.
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
% Extract EC tokens of the form a.b.c.X where X can be digits or letters.
% We then normalize (old->new mapping + keep numeric only) in normalizeEcTokens().
    s = char(string(s));
    m = regexp(s, '(?<!\d)(\d+)\.(\d+)\.(\d+)\.([A-Za-z0-9]+)(?!\d)', 'match');
    if isempty(m)
        tokens = strings(0, 1);
    else
        tokens = unique(string(m(:)));
    end
end

function ecs = normalizeEcTokens(tokens, ecMap)
% Normalize EC tokens:
%   1) apply old->new mapping if provided
%   2) keep only numeric 4-level ECs (digits.digits.digits.digits)
    tokens = string(tokens(:));
    out = strings(0, 1);
    for i = 1:numel(tokens)
        t = normalizeText(tokens(i));
        if t == ""
            continue;
        end
        if ~isempty(ecMap) && isa(ecMap, 'containers.Map') && isKey(ecMap, char(t))
            mapped = string(ecMap(char(t)));
            out = [out; mapped(:)]; %#ok<AGROW>
        else
            out = [out; t]; %#ok<AGROW>
        end
    end
    out = unique(out);
    keep = ~cellfun(@isempty, regexp(cellstr(out), '^\d+\.\d+\.\d+\.\d+$', 'once'));
    ecs = unique(out(keep));
end

function ecMap = loadEcOldNewMap(mapPath)
% Load metaMet EC old->new mapping into a containers.Map where each key maps
% to a string array (since one old EC can map to multiple new ECs).
    T = readtable(mapPath, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    v = T.Properties.VariableNames;
    oldIdx = find(strcmpi(v, 'old ec number'), 1);
    newIdx = find(strcmpi(v, 'new number'), 1);
    if isempty(oldIdx) || isempty(newIdx)
        % Fallback: assume first two columns
        oldIdx = 1;
        newIdx = 2;
    end

    olds = string(T{:, oldIdx});
    news = string(T{:, newIdx});

    olds = arrayfun(@normalizeText, olds);
    news = arrayfun(@normalizeText, news);
    keep = olds ~= "" & news ~= "";
    olds = olds(keep);
    news = news(keep);

    ecMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(olds)
        k = char(olds(i));
        vNew = news(i);
        if isKey(ecMap, k)
            ecMap(k) = unique([string(ecMap(k)); vNew]);
        else
            ecMap(k) = vNew;
        end
    end
end

function model = ensureNestedSubsystems(model)
% RAVEN's checkModelStruct (used by writeYAMLmodel) expects subSystems to be:
%   cell array (nRxns x 1) where each entry is a cell array of strings.
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

    out = repmat({{}}, nRxns, 1);
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
            entry = {};
        end
        entry = entry(:)';
        entry = entry(~cellfun(@(x) isempty(strtrim(char(string(x)))), entry));
        if isempty(entry)
            entry = {''};
        end
        out{i} = entry;
    end
    model.subSystems = out;
end
