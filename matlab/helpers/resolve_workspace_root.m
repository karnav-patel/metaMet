function workspaceRoot = resolve_workspace_root()
% Resolve the repository workspace root for ModelYML pipeline scripts.

workspaceRoot = getenv('MODELYML_WORKSPACE_ROOT');
if ~isempty(workspaceRoot) && isfolder(workspaceRoot)
    return;
end

workspaceRoot = pwd;
end