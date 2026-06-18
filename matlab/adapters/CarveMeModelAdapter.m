classdef CarveMeModelAdapter < ModelAdapter
    % Minimal adapter for building a GECKO ecModel from a CarveMe-derived GEM.
    % This adapter is intentionally lightweight: it points GECKO to a local
    % data folder (gecko_adapter/data) where a stub uniprot.tsv can be placed.
    methods
        function obj = CarveMeModelAdapter()
            workspaceRoot = resolve_workspace_root();
            layout = resolve_modelyml_layout(workspaceRoot);
            obj.params.path = layout.lightAdapterDir;
            obj.params.convGEM = layout.draftModelXml;

            obj.params.sigma = 0.5;
            obj.params.Ptot = 0.5;
            obj.params.f = 0.5;
            obj.params.gR_exp = 0.1;
            obj.params.org_name = 'CarveMe model';

            obj.params.complex.taxonomicID = [];

            % Leave external DB IDs empty to avoid accidental downloads.
            obj.params.kegg.ID = '';
            obj.params.kegg.geneID = 'kegg';

            obj.params.uniprot.type = 'taxonomy';
            obj.params.uniprot.ID = '';
            obj.params.uniprot.geneIDfield = 'Gene';
            obj.params.uniprot.reviewed = false;

            obj.params.c_source = '';
            obj.params.bioRxn = '';

            % Must match one of model.compNames in your GEM (RAVEN import keeps names).
            obj.params.enzyme_comp = 'cytosol';
        end

        function [spont,spontRxnNames] = getSpontaneousReactions(obj, model)
            if isfield(model, 'rxnNames')
                spont = contains(lower(model.rxnNames), 'spontaneous');
                spontRxnNames = model.rxnNames(spont);
            else
                spont = false(numel(model.rxns), 1);
                spontRxnNames = {};
            end
        end
    end
end
