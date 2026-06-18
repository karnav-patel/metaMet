% SBML -> YAML roundtrip using RAVEN or COBRA Toolbox
% - Input:  model.xml (SBML)
% - Output: model.yml (RAVEN-style YAML)
%
% Prereqs:
%   - RAVEN Toolbox on MATLAB path (preferred) OR COBRA Toolbox
%   - For YAML I/O: RAVEN functions writeYAMLmodel/readYAMLmodel

clearvars;

workspaceRoot = resolve_workspace_root();
layout = resolve_modelyml_layout(workspaceRoot);

sbmlPath = layout.draftModelXml;
yamlPath = layout.draftModelYaml;

if ~isfolder(layout.draftsDir)
    mkdir(layout.draftsDir);
end

if ~isfile(sbmlPath)
    error('SBML file not found: %s', sbmlPath);
end

% 1) Load SBML into a COBRA/RAVEN model struct
if exist('importModel', 'file') == 2
    % RAVEN import
    model = importModel(sbmlPath);
elseif exist('readCbModel', 'file') == 2
    % COBRA import fallback
    model = readCbModel(sbmlPath);
else
    error(['Neither RAVEN importModel nor COBRA readCbModel is on the MATLAB path.\n' ...
           'Add RAVEN (preferred) or COBRA Toolbox to the path, then re-run.']);
end

% 2) Ensure required fields are correctly typed
% NOTE: In many RAVEN versions, subSystems is expected to be a *cell array of
% cell arrays* (one entry per reaction; each entry can hold multiple
% subsystems). A plain cellstr will fail checkModelStruct.
if ~isfield(model, 'subSystems') || isempty(model.subSystems)
    if isfield(model, 'rxns')
        model.subSystems = repmat({{''}}, numel(model.rxns), 1);
    else
        model.subSystems = {{''}};
    end
else
    if ~iscell(model.subSystems)
        tmp = cellstr(model.subSystems);
        model.subSystems = cellfun(@(s) {s}, tmp(:), 'UniformOutput', false);
    else
        % Convert cellstr -> cell-of-cells, and normalize any other entries
        model.subSystems = model.subSystems(:);
        model.subSystems = cellfun(@(x) normalizeSubsystemEntry(x), model.subSystems, 'UniformOutput', false);
    end
end

% RAVEN checkModelStruct is strict about duplicate metabolite *names* within the
% same compartment. This is common in SBMLs (names are descriptive, IDs are
% unique) but it blocks YAML export. Disambiguate names for duplicates by
% appending the unique metabolite ID.
if isfield(model, 'metNames') && isfield(model, 'metComps') && isfield(model, 'comps') && isfield(model, 'mets')
    metInComp = strcat(model.metNames, '[', model.comps(model.metComps), ']');
    [uVals, ~, uIdx] = unique(metInComp);
    counts = accumarray(uIdx, 1);
    dupVals = uVals(counts > 1);
    if ~isempty(dupVals)
        dupMask = ismember(metInComp, dupVals);
        dupIdxs = find(dupMask);
        for j = 1:numel(dupIdxs)
            k = dupIdxs(j);
            model.metNames{k} = sprintf('%s (%s)', model.metNames{k}, model.mets{k});
        end
    end
end

% 3) Export RAVEN-style YAML
if exist('writeYAMLmodel', 'file') ~= 2
    error('writeYAMLmodel not found on path. Install/enable RAVEN YAML I/O.');
end
writeYAMLmodel(model, yamlPath);

% 4) Read back (sanity check)
if exist('readYAMLmodel', 'file') ~= 2
    error('readYAMLmodel not found on path. Install/enable RAVEN YAML I/O.');
end
model2 = readYAMLmodel(yamlPath);

fprintf('Wrote: %s\n', yamlPath);
if isfield(model, 'rxns')
    fprintf('Original: %d rxns, %d mets\n', numel(model.rxns), numel(model.mets));
end
if isfield(model2, 'rxns')
    fprintf('Reloaded: %d rxns, %d mets\n', numel(model2.rxns), numel(model2.mets));
end

function entry = normalizeSubsystemEntry(x)
% Ensure one reaction's subSystems entry is a cell array of char vectors.
    if isempty(x)
        entry = {''};
        return;
    end
    if iscell(x)
        % If it's a nested cell already, convert any strings/chars inside to char
        entry = x(:)';
        entry = cellfun(@(y) char(string(y)), entry, 'UniformOutput', false);
        return;
    end
    if isstring(x)
        entry = cellstr(x);
        return;
    end
    if ischar(x)
        entry = {x};
        return;
    end
    % Fallback: stringify unknown types
    entry = {char(string(x))};
end
