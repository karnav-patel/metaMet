% Quick QC for GECKO/RAVEN YAML models in this workspace.
% - Checks that each model can solve an LP (FBA)
% - Reports objective value and basic ec/kcat coverage (if present)

clearvars;

here = resolve_workspace_root();
layout = resolve_modelyml_layout(here);

% Add GECKO (solveLP lives in RAVEN/GECKO toolchain in this workspace)
geckoRoot = resolve_gecko_root(here);
if isfolder(geckoRoot)
    addpath(genpath(fullfile(geckoRoot, 'src')));
end

if exist('readYAMLmodel', 'file') ~= 2
    error('readYAMLmodel not found on MATLAB path. Ensure RAVEN is available (GECKO bundles required functions under GECKO-main/src).');
end
if exist('solveLP', 'file') ~= 2
    error('solveLP not found on MATLAB path.');
end

checks = {
    struct('label', 'model.xml',      'path', layout.draftModelXml,      'loader', 'import')
    struct('label', 'ecModel.yml',    'path', layout.lightEcModel,    'loader', 'yaml')
    struct('label', 'ecModel_light.yml', 'path', layout.lightEcModelLight, 'loader', 'yaml')
    struct('label', 'ecModel_kcat.yml', 'path', layout.lightEcModelKcat, 'loader', 'yaml')
};

for i = 1:numel(checks)
    check = checks{i};
    if ~isfile(check.path)
        fprintf('%s: missing\n', check.label);
        continue;
    end

    fprintf('\n== %s ==\n', check.label);
    [m, loadNote] = loadModelForQC(check);
    if ~isempty(loadNote)
        fprintf('%s\n', loadNote);
    end

    if isfield(m, 'c')
        fprintf('Objective nnz(c): %d\n', nnz(m.c));
    end

    if any(strcmp(m.rxns, 'prot_pool_exchange'))
        fprintf('prot_pool_exchange lb: %g\n', m.lb(strcmp(m.rxns, 'prot_pool_exchange')));
    end
    fprintf('usage_prot reactions: %d\n', sum(startsWith(string(m.rxns), 'usage_prot_')));

    if isfield(m, 'ec') && isfield(m.ec, 'kcat')
        k = m.ec.kcat;
        fprintf('ec.rxns: %d\n', numel(k));
        fprintf('kcat_nonzero: %d\n', sum(k > 0));
    end

    [isValid, validationMsg] = validateModelForLP(m);
    if ~isValid
        fprintf('FBA skipped: %s\n', validationMsg);
        continue;
    end

    try
        sol = solveLP(m);
        fprintf('FBA stat=%d obj=%g\n', sol.stat, sol.f);
    catch ME
        fprintf('FBA failed: %s\n', ME.message);
    end
end

function [model, note] = loadModelForQC(check)
    note = '';
    switch check.loader
        case 'import'
            if exist('importModel', 'file') ~= 2
                error('importModel not found on MATLAB path, so %s could not be loaded.', check.label);
            end
            model = importModel(check.path);
            note = 'Loaded via importModel to avoid draft-YAML round-trip inconsistencies.';
        case 'yaml'
            model = readYAMLmodel(check.path);
        otherwise
            error('Unsupported QC loader: %s', check.loader);
    end
end

function [tf, msg] = validateModelForLP(model)
    tf = true;
    msg = '';

    if ~isfield(model, 'S') || ~isfield(model, 'rxns') || ~isfield(model, 'mets')
        tf = false;
        msg = 'missing one or more required fields: S, rxns, mets';
        return;
    end
    if ~isreal(model.S)
        tf = false;
        msg = 'stoichiometric matrix S is not real-valued';
        return;
    end
    if size(model.S, 1) ~= numel(model.mets)
        tf = false;
        msg = sprintf('S has %d rows but model lists %d metabolites', size(model.S, 1), numel(model.mets));
        return;
    end
    if size(model.S, 2) ~= numel(model.rxns)
        tf = false;
        msg = sprintf('S has %d columns but model lists %d reactions', size(model.S, 2), numel(model.rxns));
        return;
    end
end
