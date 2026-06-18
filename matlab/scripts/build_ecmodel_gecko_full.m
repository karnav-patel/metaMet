% Build a full GECKO ecModel from the CarveMe draft and prepare the
% additional metadata required for full-model simulations.
%
% Outputs:
%   - ecModel_full.yml
%   - gecko_full_adapter/data/ComplexPortal.json (if downloadable)
%   - gecko_full_adapter/data/foundComplex.tsv (if complex data matched)
%   - gecko_full_adapter/data/proposedComplex.tsv (if partial matches exist)

clearvars;

workspaceRoot = resolve_workspace_root();
geckoRoot = resolve_gecko_root(workspaceRoot);
layout = resolve_modelyml_layout(workspaceRoot);

addpath(genpath(fullfile(geckoRoot, 'src')));
addpath(workspaceRoot);

if exist('readYAMLmodel', 'file') ~= 2
    error('readYAMLmodel not found on path. Ensure RAVEN Toolbox is installed and on path.');
end

if ~isfile(layout.draftModelYaml) && ~isfile(layout.draftModelXml)
    error('Missing both %s and %s. Run the CarveMe + SBML->YAML steps first.', layout.draftModelXml, layout.draftModelYaml);
end

if exist('prepare_gecko_full_adapter', 'file') == 2
    prepare_gecko_full_adapter;
end

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

model = ensureNestedSubsystems(model);

if exist('ModelAdapterManager', 'class') ~= 8
    error('ModelAdapterManager class not found. Verify GECKO is on the MATLAB path.');
end

adapterPath = fullfile(layout.matlabAdaptersDir, 'CarveMeFullModelAdapter.m');
if ~isfile(adapterPath)
    error('Missing full-model adapter file: %s', adapterPath);
end
modelAdapter = ModelAdapterManager.setDefault(adapterPath, true);
params = modelAdapter.getParameters();

if exist('makeEcModel', 'file') ~= 2
    error('makeEcModel not found on path. Verify GECKO-main/src is on the MATLAB path.');
end

[ecModel, noUniprot] = makeEcModel(model, false, modelAdapter);
fprintf('Full ecModel build completed. Genes without UniProt match: %d\n', numel(noUniprot));

if exist('applyComplexData', 'file') == 2
    try
        [ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel, [], modelAdapter, true);
        if ~isempty(foundComplex)
            writetable(foundComplex, fullfile(params.path, 'data', 'foundComplex.tsv'), 'FileType', 'text', 'Delimiter', '\t');
        end
        if ~isempty(proposedComplex)
            writetable(proposedComplex, fullfile(params.path, 'data', 'proposedComplex.tsv'), 'FileType', 'text', 'Delimiter', '\t');
        end
    catch ME
        warning('Complex data integration was skipped: %s', ME.message);
    end
end

ecModel = ensureNestedSubsystems(ecModel);
if ~isfolder(layout.fullModelsDir)
    mkdir(layout.fullModelsDir);
end
outYaml = layout.fullEcModel;
if exist('writeYAMLmodel', 'file') == 2
    writeYAMLmodel(ecModel, outYaml);
else
    error('writeYAMLmodel not found on path.');
end

fprintf('Wrote %s\n', outYaml);
syncFileToTutorialOutput(outYaml, 'ecModel_full.yml');
if isfield(ecModel, 'ec')
    disp('full ecModel.ec exists.');
else
    warning('full ecModel.ec does not exist after makeEcModel; check GECKO version and inputs.');
end

function syncFileToTutorialOutput(sourcePath, targetName)
    outDir = fullfile(resolve_gecko_root(resolve_workspace_root()), 'tutorials', 'full_ecModel', 'output');
    if ~isfolder(outDir)
        return;
    end
    copyfile(sourcePath, fullfile(outDir, targetName));
    fprintf('Synced %s into %s\n', targetName, outDir);
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