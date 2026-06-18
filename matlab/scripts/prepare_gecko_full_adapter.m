% Prepare a GECKO full-model adapter folder with the local files required
% to run a full ecModel workflow for this CarveMe-derived organism.
%
% Outputs in gecko_full_adapter/data:
%   - uniprot.tsv            (Entry/Gene/EC/Mass/Sequence)
%   - uniprotConversion.tsv  (model gene -> UniProt accession)
%   - pseudoRxns.tsv         (non-enzymatic reactions to ignore in standard kcat)

clearvars;

workspaceRoot = resolve_workspace_root();
layout = resolve_modelyml_layout(workspaceRoot);

adapterDir = layout.fullAdapterDir;
dataDir = layout.fullAdapterDataDir;
uniprotTsv = fullfile(dataDir, 'uniprot.tsv');
conversionTsv = fullfile(dataDir, 'uniprotConversion.tsv');
pseudoRxnsTsv = fullfile(dataDir, 'pseudoRxns.tsv');
manifestTxt = fullfile(dataDir, 'required_inputs_summary.txt');

geckoRoot = resolve_gecko_root(workspaceRoot);
addpath(genpath(fullfile(geckoRoot, 'src')));

faaCandidates = {
    layout.genomeFaa
    layout.sanitizedFasta
};

faaPath = '';
for i = 1:numel(faaCandidates)
    if isfile(faaCandidates{i})
        faaPath = faaCandidates{i};
        break;
    end
end
if isempty(faaPath)
    error('Could not find a protein FASTA at genome_for_carveme.faa or genome.faa.');
end

if exist('readYAMLmodel', 'file') == 2 && isfile(layout.draftModelYaml)
    model = readYAMLmodel(layout.draftModelYaml);
elseif exist('importModel', 'file') == 2 && isfile(layout.draftModelXml)
    model = importModel(layout.draftModelXml);
else
    error('Need model.yml + readYAMLmodel OR model.xml + importModel.');
end

if ~isfield(model, 'genes') || isempty(model.genes)
    error('Model has no genes; GECKO full-model reconstruction requires gene-associated reactions.');
end
if exist('calculateMW', 'file') ~= 2
    error('GECKO calculateMW not found on path. Add GECKO-main/src to MATLAB path first.');
end

if ~isfolder(dataDir)
    mkdir(dataDir);
end

seqMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
upMap = containers.Map('KeyType', 'char', 'ValueType', 'char');
ecMap = containers.Map('KeyType', 'char', 'ValueType', 'char');

fidFaa = fopen(faaPath, 'rt');
if fidFaa < 0
    error('Could not read %s', faaPath);
end

currentGene = '';
currentUp = '';
currentEc = '';
currentSeq = '';
while true
    ln = fgetl(fidFaa);
    if ~ischar(ln)
        break;
    end
    if startsWith(ln, '>')
        if ~isempty(currentGene)
            seqMap(currentGene) = currentSeq;
            upMap(currentGene) = currentUp;
            ecMap(currentGene) = currentEc;
        end
        headerToken = strtrim(ln(2:end));
        sp = strfind(headerToken, ' ');
        if ~isempty(sp)
            headerToken = headerToken(1:sp(1)-1);
        end
        currentGene = normalizeGeneId(extractGeneToken(headerToken));
        currentUp = extractHeaderField(headerToken, 'UP');
        currentEc = extractHeaderField(headerToken, 'EC');
        currentSeq = '';
    else
        currentSeq = [currentSeq strtrim(ln)]; %#ok<AGROW>
    end
end
if ~isempty(currentGene)
    seqMap(currentGene) = currentSeq;
    upMap(currentGene) = currentUp;
    ecMap(currentGene) = currentEc;
end
fclose(fidFaa);

genes = model.genes(:);
found = 0;
seenEntries = containers.Map('KeyType', 'char', 'ValueType', 'double');

fid = fopen(uniprotTsv, 'wt');
if fid < 0
    error('Could not write %s', uniprotTsv);
end
fprintf(fid, 'Entry\tGene\tEC\tMass\tSequence\n');

fidConv = fopen(conversionTsv, 'wt');
if fidConv < 0
    fclose(fid);
    error('Could not write %s', conversionTsv);
end
fprintf(fidConv, 'model_gene\tuniprot_id\n');

for i = 1:numel(genes)
    g = genes{i};
    if isempty(g)
        continue;
    end
    seq = '';
    mw = 0;
    upId = ['FASTA_' g];
    ecVal = '';
    if isKey(seqMap, g)
        seq = seqMap(g);
        ecVal = ecMap(g);
        if isKey(upMap, g) && ~isempty(upMap(g))
            upId = upMap(g);
        end
        if ~isempty(seq)
            mw = round(calculateMW(seq));
            found = found + 1;
        end
    end
    upId = uniquifyEntry(upId, g, seenEntries);
    fprintf(fid, '%s\t%s\t%s\t%d\t%s\n', upId, g, ecVal, mw, seq);
    fprintf(fidConv, '%s\t%s\n', g, upId);
end

fclose(fid);
fclose(fidConv);

writePseudoRxns(model, pseudoRxnsTsv);

fidManifest = fopen(manifestTxt, 'wt');
if fidManifest >= 0
    fprintf(fidManifest, 'Generated inputs for full GECKO workflow\n');
    fprintf(fidManifest, 'Adapter folder: %s\n', adapterDir);
    fprintf(fidManifest, 'Input FASTA: %s\n', faaPath);
    fprintf(fidManifest, 'uniprot.tsv: local UniProt-like database with real UP accessions where available\n');
    fprintf(fidManifest, 'uniprotConversion.tsv: model gene -> UniProt accession mapping\n');
    fprintf(fidManifest, 'pseudoRxns.tsv: objective / pseudo reactions excluded from standard kcat assignment\n');
    fprintf(fidManifest, 'Optional runtime caches created later if available: ComplexPortal.json, kegg.tsv\n');
    fclose(fidManifest);
end

fprintf('Wrote full-model UniProt-like DB: %s\n', uniprotTsv);
fprintf('Wrote gene->UniProt mapping: %s\n', conversionTsv);
fprintf('Wrote pseudo reaction list: %s\n', pseudoRxnsTsv);
fprintf('Matched sequences for %d/%d model genes (from %s).\n', found, numel(genes), faaPath);

function fieldVal = extractHeaderField(headerToken, fieldName)
    expr = [fieldName '=([^|]+)'];
    token = regexp(headerToken, expr, 'tokens', 'once');
    if isempty(token)
        fieldVal = '';
    else
        fieldVal = strtrim(token{1});
    end
end

function geneToken = extractGeneToken(headerToken)
    parts = strsplit(headerToken, '|');
    geneToken = strtrim(parts{1});
end

function out = normalizeGeneId(headerToken)
    out = headerToken;
    out = strrep(out, '|EC=', '_EC_');
    out = strrep(out, '|', '_');
    out = strrep(out, 'EC=', 'EC_');
    out = regexprep(out, '[^A-Za-z0-9_]', '_');
    out = regexprep(out, '_+', '_');
    out = regexprep(out, '^_+', '');
    out = regexprep(out, '_+$', '');
end

function writePseudoRxns(model, outPath)
    pseudo = strings(0,1);
    if isfield(model, 'c')
        pseudo = [pseudo; string(model.rxns(model.c ~= 0))]; %#ok<AGROW>
    end
    if isfield(model, 'rxnNames')
        idx = contains(lower(string(model.rxnNames)), 'growth') | contains(lower(string(model.rxnNames)), 'biomass');
        pseudo = [pseudo; string(model.rxns(idx))]; %#ok<AGROW>
    end
    pseudo = unique(pseudo(pseudo ~= ""));

    fid = fopen(outPath, 'wt');
    if fid < 0
        error('Could not write %s', outPath);
    end
    for i = 1:numel(pseudo)
        fprintf(fid, '%s\tpseudo_or_objective\n', pseudo(i));
    end
    fclose(fid);
end

function entryId = uniquifyEntry(entryId, geneId, seenEntries)
    if ~isKey(seenEntries, entryId)
        seenEntries(entryId) = 1;
        return;
    end
    seenEntries(entryId) = seenEntries(entryId) + 1;
    entryId = sprintf('%s__%s', entryId, geneId);
end