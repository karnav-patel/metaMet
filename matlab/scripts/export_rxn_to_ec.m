% Export reaction -> EC mapping to CSV
% Reads SBML or RAVEN YAML and produces a normalized table with one EC per row.

clearvars;

workspaceRoot = resolve_workspace_root();
layout = resolve_modelyml_layout(workspaceRoot);

% Prefer YAML if present (already normalized for RAVEN), otherwise SBML.
yamlPath = layout.draftModelYaml;
sbmlPath = layout.draftModelXml;
outPath  = layout.rxnToEcCsv;

if ~isfolder(layout.draftsDir)
    mkdir(layout.draftsDir);
end

if exist('readYAMLmodel', 'file') == 2 && isfile(yamlPath)
    model = readYAMLmodel(yamlPath);
elseif exist('importModel', 'file') == 2 && isfile(sbmlPath)
    model = importModel(sbmlPath);
else
    error('Could not load model: need model.yml + readYAMLmodel OR model.xml + importModel.');
end

if ~isfield(model, 'eccodes')
    error('Model does not have an eccodes field; cannot export rxn->EC mapping.');
end

rxnIds = model.rxns(:);
ecRaw  = model.eccodes(:);

rowsRxn = {};
rowsEc  = {};
rowsEcRaw = {};

for i = 1:numel(rxnIds)
    raw = ecRaw{i};
    if isempty(raw)
        continue;
    end
    rawStr = char(string(raw));
    parts = split(string(rawStr), ';');
    parts = strtrim(parts);
    parts(parts == "") = [];
    for j = 1:numel(parts)
        rowsRxn(end+1,1) = rxnIds(i);
        rowsEc(end+1,1) = {char(parts(j))};
        rowsEcRaw(end+1,1) = {rawStr};
    end
end

T = table(rowsRxn, rowsEc, rowsEcRaw, 'VariableNames', {'rxn_id','ec_number','ec_raw'});
T = unique(T);

writetable(T, outPath);
fprintf('Wrote %s (%d rows)\n', outPath, height(T));
