classdef CarveMeFullModelAdapter < ModelAdapter
    % Model adapter for a full GECKO reconstruction from a CarveMe draft.
    % Uses a dedicated gecko_full_adapter/ folder with local caches and
    % model-specific metadata generated from the FASTA headers.
    methods
        function obj = CarveMeFullModelAdapter()
            workspaceRoot = resolve_workspace_root();
            layout = resolve_modelyml_layout(workspaceRoot);
            obj.params.path = layout.fullAdapterDir;
            obj.params.convGEM = layout.draftModelXml;

            obj.params.sigma = 0.5;
            obj.params.Ptot = 0.5;
            obj.params.f = 0.5;
            obj.params.gR_exp = 0.1;
            obj.params.org_name = 'Escherichia coli (strain K-12)';

            obj.params.complex.taxonomicID = 83333;

            obj.params.kegg.ID = 'eco';
            obj.params.kegg.geneID = 'UniProt';

            obj.params.uniprot.type = 'taxonomy';
            obj.params.uniprot.ID = '83333';
            obj.params.uniprot.geneIDfield = 'Gene';
            obj.params.uniprot.reviewed = false;

            obj.params.c_source = 'R_EX_glc__D_e';
            obj.params.bioRxn = 'Growth';
            obj.params.enzyme_comp = 'cytosol';
        end

        function genes = getUniprotCompatibleGenes(~, inGenes)
            genes = string(inGenes(:));
        end

        function uniprotIDs = getUniprotIDsFromTable(obj, modelGenes)
            uniprotIDs = obj.mapGenesToUniprot(modelGenes);
        end

        function [spont, spontRxnNames] = getSpontaneousReactions(obj, model)
            spontFile = fullfile(obj.params.path, 'data', 'spontaneousReactions.tsv');
            if isfile(spontFile)
                fID = fopen(spontFile, 'r');
                fileData = textscan(fID, '%q %q', 'Delimiter', '\t', 'HeaderLines', 1);
                fclose(fID);
                allRxns = fileData{1};
                spontMask = logical(str2double(fileData{2}));
                [spont, ~] = ismember(model.rxns, allRxns(spontMask));
                spontRxnNames = model.rxns(spont);
                return;
            end

            if isfield(model, 'rxnNames')
                spont = contains(lower(model.rxnNames), 'spontaneous');
                spontRxnNames = model.rxnNames(spont);
            else
                spont = false(numel(model.rxns), 1);
                spontRxnNames = {};
            end
        end
    end

    methods (Access = private)
        function mapped = mapGenesToUniprot(obj, genes)
            genes = cellstr(string(genes(:)));
            mapped = string(genes);
            conversionPath = fullfile(obj.params.path, 'data', 'uniprotConversion.tsv');
            if ~isfile(conversionPath)
                return;
            end

            fID = fopen(conversionPath, 'r');
            data = textscan(fID, '%q %q', 'Delimiter', '\t', 'HeaderLines', 1);
            fclose(fID);
            if isempty(data{1})
                return;
            end

            keyVals = data{1};
            mapVals = data{2};
            [found, idx] = ismember(genes, keyVals);
            mapped(found) = string(mapVals(idx(found)));
        end
    end
end