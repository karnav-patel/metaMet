function layout = resolve_modelyml_layout(workspaceRoot)
% Resolve the organized ModelYML workspace layout used by the Python driver.

if nargin < 1 || isempty(workspaceRoot)
    workspaceRoot = resolve_workspace_root();
end

inputsDir = getenv('MODELYML_INPUTS_DIR');
if isempty(inputsDir)
    inputsDir = fullfile(workspaceRoot, 'inputs');
end

artifactsDir = getenv('MODELYML_ARTIFACTS_DIR');
if isempty(artifactsDir)
    artifactsDir = fullfile(workspaceRoot, 'artifacts');
end

layout = struct();
layout.workspaceRoot = workspaceRoot;
layout.inputsDir = inputsDir;
layout.artifactsDir = artifactsDir;
layout.draftsDir = fullfile(artifactsDir, 'drafts');
layout.logsDir = fullfile(artifactsDir, 'logs');
layout.reportsDir = fullfile(artifactsDir, 'reports');
layout.collectionsDir = fullfile(artifactsDir, 'collections');
layout.modelsDir = fullfile(artifactsDir, 'models');
layout.lightModelsDir = fullfile(layout.modelsDir, 'light');
layout.fullModelsDir = fullfile(layout.modelsDir, 'full');
layout.adaptersDir = fullfile(artifactsDir, 'adapters');
layout.lightAdapterDir = fullfile(layout.adaptersDir, 'light');
layout.fullAdapterDir = fullfile(layout.adaptersDir, 'full');
layout.lightAdapterDataDir = fullfile(layout.lightAdapterDir, 'data');
layout.fullAdapterDataDir = fullfile(layout.fullAdapterDir, 'data');
layout.metaMetModelsDir = fullfile(workspaceRoot, 'metaMet', 'data_modeling', 'models');
layout.genomeFaa = fullfile(inputsDir, 'genome.faa');
layout.exampleGenomeFaa = fullfile(inputsDir, 'genome_ABC.faa');
layout.sanitizedFasta = fullfile(layout.draftsDir, 'genome_for_carveme.faa');
layout.sanitizedFastaMap = fullfile(layout.draftsDir, 'genome_for_carveme_header_map.tsv');
layout.draftModelXml = fullfile(layout.draftsDir, 'model.xml');
layout.draftModelYaml = fullfile(layout.draftsDir, 'model.yml');
layout.nonGapfillModelXml = fullfile(layout.draftsDir, 'model_nogapfill.xml');
layout.rxnToEcCsv = fullfile(layout.draftsDir, 'rxn_to_ec.csv');
layout.lightEcModel = fullfile(layout.lightModelsDir, 'ecModel.yml');
layout.lightEcModelLight = fullfile(layout.lightModelsDir, 'ecModel_light.yml');
layout.lightEcModelKcat = fullfile(layout.lightModelsDir, 'ecModel_kcat.yml');
layout.fullEcModel = fullfile(layout.fullModelsDir, 'ecModel_full.yml');
layout.fullEcModelKcat = fullfile(layout.fullModelsDir, 'ecModel_full_kcat.yml');
layout.matlabRoot = fullfile(workspaceRoot, 'matlab');
layout.matlabScriptsDir = fullfile(layout.matlabRoot, 'scripts');
layout.matlabHelpersDir = fullfile(layout.matlabRoot, 'helpers');
layout.matlabAdaptersDir = fullfile(layout.matlabRoot, 'adapters');
end
