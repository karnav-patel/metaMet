% Build a GECKO ecModel from a RAVEN GEM YAML
% Inputs:
%   - model.yml (RAVEN YAML exported previously)
% Outputs:
%   - ecModel.yml (full GECKO/RAVEN YAML containing .ec structure)
%   - ecModel_light.yml (light GECKO/RAVEN YAML containing .ec structure)

clearvars;

workspaceRoot = resolve_workspace_root();
geckoRoot = resolve_gecko_root(workspaceRoot);
layout = resolve_modelyml_layout(workspaceRoot);

% Add GECKO to path (and its dependencies if they live inside)
addpath(genpath(fullfile(geckoRoot, 'src')));

% Ensure our local adapter class is visible
addpath(workspaceRoot);

if exist('readYAMLmodel', 'file') ~= 2
    error('readYAMLmodel not found on path. Ensure RAVEN Toolbox is installed and on path.');
end

if ~isfile(layout.draftModelYaml) && ~isfile(layout.draftModelXml)
    error('Missing both %s and %s. Run the CarveMe + SBML->YAML steps first.', layout.draftModelXml, layout.draftModelYaml);
end

% Prepare adapter folder + stub UniProt DB (avoids downloads for the first run)
if exist('prepare_gecko_adapter', 'file') == 2
    prepare_gecko_adapter;
end

% Re-load the source model after prepare_gecko_adapter, because that script
% runs in the caller workspace and populates its own `model` variable.
% Prefer direct SBML import when available. Some gap-filled CarveMe YAML
% round-trips can retain reaction metadata but reload an S matrix with fewer
% columns than reactions, which breaks GECKO makeEcModel.
gemYaml = layout.draftModelYaml;
sbmlPath = layout.draftModelXml;
if exist('importModel', 'file') == 2 && isfile(sbmlPath)
    model = importModel(sbmlPath);
elseif isfile(gemYaml)
    model = readYAMLmodel(gemYaml);
else
    error('Could not load a source model. Need importModel + model.xml or readYAMLmodel + model.yml.');
end

if isfield(model, 'S') && isfield(model, 'rxns') && size(model.S, 2) ~= numel(model.rxns)
    error('Model import is inconsistent: S has %d columns but model has %d reactions.', size(model.S, 2), numel(model.rxns));
end

% GECKO/RAVEN YAML I/O expects subSystems to be:
%   cell array (nRxns x 1) where each entry is a cell array of strings.
% This is a common failure point when models were imported from SBML or
% produced by other toolchains.
model = ensureNestedSubsystems(model);

% Register default adapter (required by makeEcModel)
if exist('ModelAdapterManager', 'class') ~= 8
    error('ModelAdapterManager class not found. Verify GECKO is on the MATLAB path.');
end

adapterPath = fullfile(layout.matlabAdaptersDir, 'CarveMeModelAdapter.m');
if ~isfile(adapterPath)
    error('Missing adapter file: %s', adapterPath);
end
modelAdapter = ModelAdapterManager.setDefault(adapterPath, true);

% Build ecModel skeleton
if exist('makeEcModel', 'file') ~= 2
    error('makeEcModel not found on path. Verify GECKO-main/src is on the MATLAB path.');
end

% makeEcModel signature may vary slightly across GECKO versions.
% Build both the full and light GECKO variants from the same base model.
baseModel = model;
fullEcModel = makeEcModel(baseModel, false, modelAdapter);
lightEcModel = makeEcModel(baseModel, true, modelAdapter);

% Save in YAML. GECKO provides saveEcModel; fall back to writeYAMLmodel if needed.
if ~isfolder(layout.lightModelsDir)
    mkdir(layout.lightModelsDir);
end
fullOutYaml = layout.lightEcModel;
lightOutYaml = layout.lightEcModelLight;
if exist('writeYAMLmodel', 'file') == 2
    writeYAMLmodel(fullEcModel, fullOutYaml);
    writeYAMLmodel(lightEcModel, lightOutYaml);
else
    error('writeYAMLmodel not found on path.');
end

% Sanity checks
fprintf('Wrote %s\n', fullOutYaml);
if isfield(fullEcModel, 'ec')
    disp('full ecModel.ec exists.');
else
    warning('full ecModel.ec does not exist after makeEcModel; check GECKO version and inputs.');
end

fprintf('Wrote %s\n', lightOutYaml);
if isfield(lightEcModel, 'ec')
    disp('light ecModel.ec exists.');
else
    warning('light ecModel.ec does not exist after makeEcModel; check GECKO version and inputs.');
end

function model = ensureNestedSubsystems(model)
% RAVEN's checkModelStruct (used by writeYAMLmodel) expects subSystems to be:
%   cell array (nRxns x 1) where each entry is a cell array of strings.
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

    % Ensure outer container is a cell array.
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
