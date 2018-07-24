

function hg38match = comparehg38VSAllen(hg38match, allenProbes)

for k=1:size(hg38match,1)
    % select probe
    probeName = hg38match.probeNames{k};
    probeID1 = hg38match.ID(k);
    % make sure that this probes exists in allen data and it was uniquely
    % reannotated
    if ~isempty(probeName) && sum(strcmp(probeName, allenProbes.probeNames(:)))~=0
        
        ind = strcmp(probeName, allenProbes.probeNames(:));
        probeID2 = allenProbes.NCBIgeneID(ind);
        % check if gene IDs match
        hg38match.allenID(k) = probeID2;
        % if entrezIDs match, label
        if probeID1==probeID2
            hg38match.compare(k) = 1; % 1 if matching
        elseif isnan(probeID2) && ~isnan(probeID1)
            hg38match.compare(k) = 2; % 2 if updated
        elseif isnan(probeID1)
            hg38match.compare(k) = 3; % NCBI didn't give a valid entrezID
        else
            hg38match.compare(k) = 0; % 0 if mismatching
        end
    else
        % if there is no such probes in allen data, make it NaN
        hg38match.allenID(k) = NaN;
    end

end
end
