% Prepare a minimal GECKO adapter folder with a local UniProt-like DB.
% This allows makeEcModel() to run without downloading UniProt/KEGG.
%
% Strategy:
%   - Parse genome.faa (protein FASTA)
%   - Normalize FASTA headers to match model gene IDs
%   - Build gecko_adapter/data/uniprot.tsv with gene->(sequence,MW)

clearvars;

workspaceRoot = resolve_workspace_root();
layout = resolve_modelyml_layout(workspaceRoot);

adapterDir = layout.lightAdapterDir;
dataDir    = layout.lightAdapterDataDir;
uniprotTsv = fullfile(dataDir, 'uniprot.tsv');
faaCandidates = {
    layout.sanitizedFasta
    layout.genomeFaa
};

faaPath = '';
for i = 1:numel(faaCandidates)
    if isfile(faaCandidates{i})
        faaPath = faaCandidates{i};
        break;
    end
end

if isempty(faaPath)
    faaPath = layout.genomeFaa;
end

if exist('readYAMLmodel', 'file') == 2 && isfile(layout.draftModelYaml)
    model = readYAMLmodel(layout.draftModelYaml);
elseif exist('importModel', 'file') == 2 && isfile(layout.draftModelXml)
    model = importModel(layout.draftModelXml);
else
    error('Need model.yml + readYAMLmodel OR model.xml + importModel.');
end

if ~isfield(model, 'genes') || isempty(model.genes)
    error('Model has no genes; GECKO makeEcModel requires gene-associated reactions.');
end

if ~isfolder(dataDir)
    mkdir(dataDir);
end

if exist('calculateMW', 'file') ~= 2
    error('GECKO calculateMW not found on path. Add GECKO-main/src to MATLAB path first.');
end

seqMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

if isfile(faaPath)
    fidFaa = fopen(faaPath, 'rt');
    if fidFaa < 0
        error('Could not read %s', faaPath);
    end
    currentId = '';
    currentSeq = '';
    while true
        ln = fgetl(fidFaa);
        if ~ischar(ln)
            break;
        end
        if startsWith(ln, '>')
            if ~isempty(currentId)
                seqMap(currentId) = currentSeq;
            end
            headerToken = strtrim(ln(2:end));
            sp = strfind(headerToken, ' ');
            if ~isempty(sp)
                headerToken = headerToken(1:sp(1)-1);
            end
            currentId = normalizeGeneId(headerToken);
            currentSeq = '';
        else
            currentSeq = [currentSeq strtrim(ln)];
        end
    end
    if ~isempty(currentId)
        seqMap(currentId) = currentSeq;
    end
    fclose(fidFaa);
else
    warning('Protein FASTA not found at %s; will generate an empty uniprot.tsv with MW=0.', faaPath);
end

% Always (re)write uniprot.tsv so it matches the current model genes.
fid = fopen(uniprotTsv, 'wt');
if fid < 0
    error('Could not write %s', uniprotTsv);
end
fprintf(fid, 'Entry\tGene\tEC\tMass\tSequence\n');

genes = model.genes(:);
found = 0;
for i = 1:numel(genes)
    g = genes{i};
    if isempty(g)
        continue;
    end
    entry = ['FASTA_' g];
    seq = '';
    mw = 0;
    if isKey(seqMap, g)
        seq = seqMap(g);
        if ~isempty(seq)
            mw = round(calculateMW(seq));
            found = found + 1;
        end
    end
    fprintf(fid, '%s\t%s\t\t%d\t%s\n', entry, g, mw, seq);
end
fclose(fid);

fprintf('Wrote UniProt-like DB: %s\n', uniprotTsv);
fprintf('Matched sequences for %d/%d model genes (from genome.faa).\n', found, numel(genes));

function out = normalizeGeneId(headerToken)
% Convert FASTA header token into the gene ID format used by importModel.
% Example: seq40|EC=1.1.1.113 -> seq40_EC_1_1_1_113
    out = headerToken;
    out = strrep(out, '|EC=', '_EC_');
    out = strrep(out, '|', '_');
    out = strrep(out, 'EC=', 'EC_');
    % Replace anything non-alphanumeric/underscore with underscore
    out = regexprep(out, '[^A-Za-z0-9_]', '_');
    out = regexprep(out, '_+', '_');
    out = regexprep(out, '^_+', '');
    out = regexprep(out, '_+$', '');
end
