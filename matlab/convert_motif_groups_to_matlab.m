function [void] = convert_motif_groups_to_matlab(MotifLibraryLocation,Input,loop_type,version,MotifLibraryPath)

    [modified_base_to_parent,modified_base_atom_list,modified_atom_to_parent,parent_atom_to_modified,modified_base_to_hydrogens,modified_base_to_hydrogen_coordinates] = zDefineModifiedNucleotides();  % load in the modified nucleotide definitions

	if ~exist([MotifLibraryLocation Input filesep 'release.json'],'file')
		% download the release from the RNA Motif Atlas website
		% https://rna.bgsu.edu/rna3dhub/motifs/release/il/3.92/json
		url = ['https://rna.bgsu.edu/rna3dhub/motifs/release/' lower(loop_type) '/' version '/json'];
		fprintf('pJAR3DMaster: downloading %s from the RNA Motif Atlas %s\n', Input, url)
		options = weboptions('Timeout', 45); % Set timeout to 45 seconds
		websave([MotifLibraryLocation Input filesep 'release.json'],url,options);
	end

	% read the json file into text
	[MotifLibraryLocation Input filesep 'release.json']
	json = fileread([MotifLibraryLocation Input filesep 'release.json']);

	% replace ||||P_1 with empty string because motif atlas keeps them but Matlab does not
	% may cause problems later if we download interactions from the RNA 3D Hub,
	% we'll have to strip them out again there
	json = regexprep(json,'(\|\|\|\|P_1)','');

	% load the json file into a Matlab data structure
	motif_list = jsondecode(json);

	% we don't know ahead of time how many distinct pdb ids there will be
	% make a set data structure that we can add to without duplication
	pdb_ids = {};
	% motif_structure is a cell array of structured variables
	motif_structure = cell(1,length(motif_list));
    GroupData(length(motif_list)) = struct();
    pdb_id_to_motif_and_loop = containers.Map();

    % get all pdb ids, and set up empty motif structure for each motif group
	for i = 1:length(motif_list)
		i
		motif_list{i}
		loop_ids = fieldnames(motif_list{i}.alignment);
		for j = 1:numel(loop_ids)
			loop_id = loop_ids{j};
			fields = split(loop_ids{j},'_');
			pdb_id = fields{2};
			pdb_ids = [pdb_ids pdb_id];
            if ~isKey(pdb_id_to_motif_and_loop,pdb_id)
                pdb_id_to_motif_and_loop(pdb_id) = [];
            end
            pdb_id_to_motif_and_loop(pdb_id) = [pdb_id_to_motif_and_loop(pdb_id); [i j]];
        end

		num_instances = length(motif_list{i}.alignment);
		num_positions = length(motif_list{i}.alignment.(loop_ids{1}));

		% this is the structure we need to produce, although not all of these fields
		% Candidates: [278×7 uint16]
		% File: [1×278 struct]
		% Query: [1×1 struct]
		% consensusEdge: [6×6 double]
		% Signature: 'cWW-cWW-cWW'

		motif_structure{i} = struct();
        motif_structure{i}.motid_id = motif_list{i}.motif_id;
        motif_structure{i}.Candidates = zeros(num_instances,num_positions+1,'uint16');
		motif_structure{i}.consensusEdge = zeros(num_positions,num_positions);
		motif_structure{i}.bp_signature = motif_list{i}.bp_signature;  % text string
		motif_structure{i}.Signature = motif_list{i}.bp_signature;  % text string
		motif_structure{i}.chainbreak = motif_list{i}.chainbreak;

		% create Truncate variable needed later
		Truncate = zeros(length(motif_structure{i}.chainbreak),1);
		if iscell(motif_structure{i}.chainbreak)
		  for j = 1:length(motif_structure{i}.chainbreak)
			% this code uses the index of the nucleotide after the break
			Truncate(j) = str2double(motif_structure{i}.chainbreak{j})+1;
		  end
		else
		  for j = 1:length(motif_structure{i}.chainbreak)
			% this code uses the index of the nucleotide after the break
			Truncate(j) = motif_structure{i}.chainbreak(j)+1;
		  end
		end
		motif_structure{i}.Truncate = Truncate;

        % some functions look for information in a Query field
        motif_structure{i}.Query = struct();
        motif_structure{i}.Query.Name = motif_list{i}.motif_id;
        motif_structure{i}.MaxDiffMat = ones(1,num_positions);
        for k = 1:length(motif_list{i}.chainbreak)
            cb = motif_list{i}.chainbreak{k};
            % pretty sure that you put Inf to indicate a chain break after position cb
            % this is used in xFASTACandidates.m
            motif_structure{i}.MaxDiffMat(cb) = Inf;
        end

        GroupData(i).MotifID = motif_list{i}.motif_id;
        % GroupData(i).Signature = cell();
        % GroupData(i).Signature{1} = motif_list{i}.signature;
        % GroupData(i).NumNT = num_positions;
        % GroupData(i).NumInstances = num_instances;
        % GroupData(i).Structured = 1;   % just set it this way, maybe that's OK
        % GroupData(i).NumStacks = 0;
        % GroupData(i).NumBPh = 0;
        % GroupData(i).NumBR = 0;
        % GroupData(i).NumBasepairs = 0;
        % GroupData(i).OwnSequence = cell(1,num_instances);
        % GroupData(i).OwnScore = zeros(1,num_instances);
        % GroupData(i).SequenceLengths = zeros(1,num_instances);
	end

	pdb_ids = unique(pdb_ids);
	fprintf('pJAR3DMaster: found %d pdb ids\n',length(pdb_ids))

	% loop over pdb_ids to only load each one once since that is slow

	for p = 1:length(pdb_ids)
		pdb_id_loaded = pdb_ids{p};
		File = zAddNTData(pdb_id_loaded);

        % remove the AA field from File
        File = rmfield(File,'AA');
        File = rmfield(File,'Het');

		% download full loops in the structure so we can load any bulges
		loop_filename = ['..\..\..\PDBFiles' filesep pdb_id_loaded '_loops.csv'];
		if ~exist(loop_filename,'file')
			url = ['https://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/' pdb_id_loaded];
			websave(loop_filename,url);
			fprintf('Saved %s in %s', url, loop_filename);
		end

		% read the full loops in the structure into a Matlab data structure
		fid = fopen(loop_filename);
		lines = textscan(fid, '%s', 'Delimiter', '\n');
		fclose(fid);
		lines = lines{1};
		loop_id_to_unit_ids = containers.Map('KeyType','char','ValueType','char');
        % loop_id_to_border = containers.Map('KeyType','char','ValueType','uint16');
		for a = 1:length(lines)
			% remove leading and trailing double quotes
			L = lines{a}(2:end-1);
			% replace ||||P_1 with empty string because motif atlas keeps them but Matlab does not
			% may cause problems later if we download interactions from the RNA 3D Hub,
			% we'll have to strip them out again there
			L = regexprep(L,'(\|\|\|\|P_1)','');
			% split the line into fields by ","
			fields = strsplit(L,'","');
			loop_id = fields{1};
			% map loop_id to comma-separated list of unit ids
			loop_id_to_unit_ids(loop_id) = fields{2};
            % loop_id_to_border(loop_id) = strsplit(fields{3},',');
		end

		% make a container to map unit ids to indices in the File structure
		unit_id_to_index = containers.Map('KeyType','char','ValueType','uint16');
		for n = 1:length(File.NT)
			unit_id_to_index(File.NT(n).ID) = n;
		end

		% loop over motif_list many times, even though inefficient
		for i = 1:length(motif_list)
			% motif_list{i}.motif_id
			loop_ids = fieldnames(motif_list{i}.alignment);

			% loop over loop_ids in the motif group
			for j = 1:numel(loop_ids)
				loop_id = loop_ids{j};
				fields = split(loop_id,'_');
				pdb_id = fields{2};

				if strcmp(pdb_id,pdb_id_loaded)
					fprintf('%s %s %s %s\n',motif_list{i}.motif_id,loop_id,pdb_id,pdb_id_loaded)

					% indices of core and bulged positions in the overall File.NT
					all_unit_ids = strsplit(loop_id_to_unit_ids(loop_id),",");
					all_indices = zeros(1,length(all_unit_ids));
                    % border = 0;
                    % sequence = '';
                    % position_to_border = loop_id_to_border(loop_id);
					for k = 1:length(all_unit_ids)
						all_indices(k) = map_unit_id_to_index(unit_id_to_index,all_unit_ids{k});
                        % fields = split(all_unit_ids{k},'|');
                        % unit = fields{4};
                        % if ismember(unit,{'A','C','G','U'})
                        %     parent = unit;
                        % elseif isKey(modified_base_to_parent,unit)
                        %     parent = modified_base_to_parent(unit);
                        % else
                        %     fprintf('Could not find parent for %s\n',unit)
                        %     crash
                        % end
                        % sequence = [sequence parent];

                        % border = border + str2num(position_to_border{k});

                        % % if border is a multiple of 2,
                        % if border % 2 == 0
                        %     sequence = [sequence '*'];
                        % end
					end

                    % GroupData(i).OwnSequence{j} = sequence;

                    % make sub-File including all nucleotides in each strand
					sf = zSubFile(File,all_indices);

					% map unit ids to indices in the sub-file
					sf_unit_id_to_index = containers.Map('KeyType','char','ValueType','uint16');
					for k = 1:length(sf.NT)
						sf_unit_id_to_index(sf.NT(k).ID) = k;
                        if sf.NT(k).Code > 4
                            C = find(sf.NT(k).Base == 'ACGU');    % code number of parent base
                            sf.NT(k).Code = C;                 % ignore that this is a modified base
                        end
					end

					% get indices of each unit id relative to the subfile to save in candidates
					% unit_ids in core positions in this loop
					unit_ids = motif_list{i}.alignment.(loop_id);
					candidate_indices = zeros(1,length(unit_ids)+1);
					for k = 1:length(unit_ids)
						candidate_indices(k) = map_unit_id_to_index(sf_unit_id_to_index,unit_ids{k});
					end
					% list column of candidate_indices points to the sub-file index
					candidate_indices(end) = j;

                    % make sure that the Edge matrix of basepairs includes all
                    % flanking basepairs, even if they might not have been annotated
                    % that way by Matlab code.
                    % Some basepairs are assumed to be there because they are AU or GC or GU
                    % even if they do not have good enough geometry to be annotated.
                    % first and last nucleotides make a cWW pair in every loop
                    sf.Edge(candidate_indices(1),candidate_indices(end-1)) = 1;
                    for k = 1:length(motif_list{i}.chainbreak)
                        cb = motif_list{i}.chainbreak{k};
                        sf.Edge(cb,cb+1) = 1;
                    end
					% add a row to the candidate list
					motif_structure{i}.Candidates(j,:) = candidate_indices;
					% add the sub-file to the motif group in the correct place
					motif_structure{i}.File(j) = sf;
				end
			end
		end
	end

    % save the motif group data
    % group_data_filename = [MotifLibraryLocation Input filesep upper(loop_type) '_GroupData.mat'];
    % save(group_data_filename,'GroupData','-v7.3');

	% loop over motif groups and save the individual files
	fprintf('Saving motif group .mat files\n')
	for i = 1:length(motif_list)
		% save the motif group
        fprintf('Saving %s\n',motif_list{i}.motif_id)
		fn = [MotifLibraryPath filesep motif_list{i}.motif_id '.mat'];
		Search = motif_structure{i};
		save(fn,'Search','-v7.3');
	end
