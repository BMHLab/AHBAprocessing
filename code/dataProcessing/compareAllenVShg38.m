% compare the other way around - how many of allen probes with entrez IDs
% are re-annotated 
% remove probes without entrez ID
function allenProbes = compareAllenVShg38(hg38match, allenProbes)

ind = isnan(allenProbes.NCBIgeneID);
allenProbes(ind,:) = []; 

numAllenGenes = size(allenProbes,1); 
for k=1:numAllenGenes
    % select probe
    probeName = allenProbes.probeNames{k};
    probeID1 = allenProbes.NCBIgeneID(k);
    % make sure that this probes exists in allen data and it was uniquely
    % reannotated
    if ~isempty(probeName) && sum(strcmp(probeName, hg38match.probeNames(:)))~=0
        
        ind = find(strcmp(probeName, hg38match.probeNames(:)));
        probeID2 = hg38match.ID(ind);
        % check if gene IDs match
        allenProbes.hg38(k) = probeID2;
        % if entrezIDs match, label
        if probeID1==probeID2
            allenProbes.compare(k) = 1;
        else
            allenProbes.compare(k) = 0;
        end
    else
        % if there is no such probes in allen data, make it NaN
        allenProbes.hg38(k) = NaN;
    end

end
end