function geckoRoot = resolve_gecko_root(workspaceRoot)
% Resolve the GECKO checkout used by ModelYML pipeline scripts.

if nargin < 1 || isempty(workspaceRoot)
    workspaceRoot = resolve_workspace_root();
end

envPath = getenv('GECKO_MAIN_DIR');
if ~isempty(envPath) && isfolder(envPath)
    geckoRoot = envPath;
    return;
end

candidates = {
    fullfile(workspaceRoot, 'external', 'GECKO-main')
    fullfile(workspaceRoot, 'GECKO-main')
};

geckoRoot = '';
for i = 1:numel(candidates)
    if isfolder(candidates{i})
        geckoRoot = candidates{i};
        return;
    end
end

error('Could not find a GECKO checkout. Set GECKO_MAIN_DIR or run the Python pipeline so it can clone GECKO into external/GECKO-main.');
end