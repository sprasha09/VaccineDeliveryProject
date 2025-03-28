function result = runGCsMain(varargin)
% Summary: Store summary statistics of the simulation status.

% INPUTS:
%   varargin: Parameters defining simulation conditions. 
%             Refer to documentation for "initializeParameters".

% OUTPUT:
%   result: Struct containing the following fields:

%   - param: Parameters struct. 
%     [Refer to documentation for initializeParameters.m]

%   - naive: 200x2210x7 array representing naive B cells.
%     - Dim1: GC
%     - Dim2: Lineage
%     - Dim3: [Lineage, Target, WT-aff, Variant-aff, Time of activation,
%              First 40 residues (decimal), Last 40 residues (decimal)].
%     * Convert decimal residues to binary for mutation state.

%   - gc: Information about GCs.
%     - gc.numbytime: 200xNx4 array (# of GC B cells with affinities >6, 7, 8, 9).
%         * N = (tmax*4+1)*2.
%         - Dim1: GC
%         - Dim2: [Dominant B cells (First N/2 columns), Sub B cells (Last N/2 columns)].
%         - Dim3: Time from 0 to tmax at 0.25-day intervals.
%     - gc.affbytime: Like gc.numbytime, but for 100, 90, 75, 50 percentiles affinities.
%     - gc.numbylineage: 200xNx2010 array (# of GC B cells per lineage).
%         * N = tmax/7+1.
%         - Dim1: GC
%         - Dim2: Time from 0 to tmax at 7-day intervals.
%         - Dim3: Lineage number.
%     - gc.finalgc: 200x3000x5 array of GC B cells at the last time point.
%         - Dim1: GC
%         - Dim2: B cells
%         - Dim3: [Lineage, Target, WT-aff, Var-aff, # of mutations].

%   - conc: Concentration information.
%     - conc.concarray: 4x3x(tmax*4+1) array for antigen & antibodies concentrations.
%         - Dim1: [Ag, IgM(natural), IgM(immune), IgG]
%         - Dim2: [Soluble, IC1, IC2 (for Ag); Empty, Target1, Target2 (for Abs)].
%         [Refer to updateConcentrations.m for details]
%     - conc.concarray_Epmask: 2x2x(tmax*4+1) for antigen with/without epitope masking.
%         - Dim1: [Soluble, IC-FDC]
%         - Dim2: [Dominant, Subdominant epitope].
%     - conc.Kaarray & conc.Kaarray_var: 3x2x(tmax*4+1) for antibody WT/Variant affinities.
%         - Dim1: [IgM(natural), IgM(immune), IgG]
%         - Dim2: [Dominant epitope, Subdominant epitope].

%   - output: GC output cell details.
%     Contains fields similar to "gc" but specific for PCs and Mem B cells.

%   - memoryCellsEGC: 9xN array of EGC-derived Memory cells.
%   - plasmaCellsEGC: 6xN array of EGC-derived PCs.

%   - dead: Dead PCs information.
%     - dead.plasmaCells & dead.plasmaCellsEGC: 6x(2*10^6) arrays for dead GC/EGC-derived PCs.
%     - dead.PCnum & dead.numPC: Scalars for the number of dead GC/EGC-derived PCs.

%--------------------------------------------------

%% Define file name for saving the data
% File name and location are determined based on the simulation parameters
addpath('helper_functions')
addpath('tcell_functions')
addpath('sim_initialization_functions') 
addpath('concentration_functions')

if nargin > 1 %Individual parameters are given, rather than the parameter struct
    param = initializeParameters(varargin{:});
    saveresult = 1; %Change to 0 if don't want to save result
else
    param = varargin{1};
    saveresult = 0;
end
rng(param.first);
fnm = getFileLocation(param);
tstart = tic;

%% Initialization
M_GC = param.M_GC; % Number of GCs;
N_GC_MAX = param.N_GC_MAX; % Max number of GC B cells
N_PC_MAX = param.N_PC_MAX; % Max number of PCs (gc B cells become PCs)
N_PB_MAX = param.N_PB_MAX; % Max number of PBs (naive B cells become PBs)
N_MEM_MAX = param.N_PC_MAX; % Max number of memory cells
n_res = param.n_res;

[naiveBcells, mutations] = getNaiveBcells(param); % Initialize naive B cells
gcBcells = zeros(M_GC, N_GC_MAX, param.gcfieldnum);
 %Lineage, target, affinity, variant affinity, nummut
gcMutations = zeros(M_GC, N_GC_MAX*n_res);
  %Each GC B cell mutation state is a string of 0s and 1s of length n_res
plasmaCells = zeros(M_GC, N_PC_MAX, param.pcfieldnum);
  %Lineage, target, affinity, variant affinity, nummut, activatedtime 
plasmaBlasts = zeros(M_GC, N_PB_MAX, param.pcfieldnum-1);
  %Lineage, target, affinity, variant affinity, activatedtime 
memoryCells = zeros(M_GC, N_MEM_MAX, param.memfieldnum);
  %Lineage, target, affinity, variant affinity, nummut, activatedtime, 
  %mutationstaet1, mutationstate2, uniqueCloneIndex
dead.plasmaCells = zeros(param.pcfieldnum, N_PC_MAX*M_GC);
dead.PCnum = 0;

numGC = zeros(M_GC,1); %This will track the last non-empty entries 
%of each row in gcBcells array; NOT equal to the number of GC B cells
numPB = zeros(M_GC,1);
numPC = zeros(M_GC,1);
numMC = zeros(M_GC,1);

result = [];
% Flatten the arrays to 2D for the convenience of operations
flatten2D = @(A) reshape(permute(A,[2,1,3]), size(A,2),[])';
naiveBcells = flatten2D(naiveBcells);
gcBcells = flatten2D(gcBcells);
plasmaCells = flatten2D(plasmaCells);
plasmaBlasts = flatten2D(plasmaBlasts);
memoryCells = flatten2D(memoryCells);

[agconc, abconc, Ka, param] = getDosingParameters(param);
Ka_var = Ka;

%% Load pre-existing memory cells and concentrations
if param.vaxnum>1 %Vax 2 and later
    %Load previous data, Initialize arrays
    previous = load(fnm{param.vaxnum-1}, 'result');
    previous = previous.result;
    flatten1D = @(A) reshape(permute(A,[3,2,1]), size(A,3),[]);
    memoryCellsEGC = zeros(param.memfieldnum, N_PC_MAX*M_GC);
    numMCEGC = 0;
    plasmaCellsEGC = zeros(param.pcfieldnum, N_PC_MAX*M_GC);
    numPCEGC = 0;

    % Memory cells
    GC_derived_mem = alignLeft(flatten1D(previous.output.finalmem), 1);
    EGC_derived_mem = alignLeft(previous.memoryCellsEGC, 1);
    existingMemory = [GC_derived_mem, EGC_derived_mem];
    [memoryToGC, mutations, memoryToEGC] =... %Get memory B cells 
        splitMemory(existingMemory, ... that will go into re-entry pool
        previous.output.memMutations, mutations, param);
    naiveBcells = [naiveBcells, memoryToGC];
    memoryCellsEGC(:,1:size(memoryToEGC,2)) = memoryToEGC;
    numMCEGC = numMCEGC + size(memoryToEGC,2);

    % Plasma cells
    existingPC = alignLeft(flatten1D(previous.output.finalpc),1);
    if ~isempty(existingPC)
        numPCEGC = size(existingPC,2);
        plasmaCellsEGC(:,1:numPCEGC) = existingPC;
    end
    dead.plasmaCellsEGC = zeros(param.pcfieldnum, N_PC_MAX*M_GC);
    dead.numPC = 0;

    % Load antigen concentration, antibody concentration and affinity
    agconc = squeeze(previous.conc.concarray(1,1:param.n_ep+1,end));
    abconc = squeeze(previous.conc.concarray(2:4,2:param.n_ep+1,end));
    Ka = squeeze(previous.conc.Kaarray(:,:,end));
    Ka_var = squeeze(previous.conc.Kaarray_var(:,:,end));
    
else % Vax 1
    % Memory cells
    memoryCellsEGC = zeros(param.memfieldnum,N_PC_MAX*M_GC);
    numMCEGC = 0;
    % Plasma cells
    plasmaCellsEGC = zeros(param.pcfieldnum,N_PC_MAX*M_GC);
    %Lineage, target, affinity, variant affinity, nummut, differentiation time 
    numPCEGC = 0;
    dead.plasmaCellsEGC = zeros(param.pcfieldnum,N_PC_MAX*M_GC);
    dead.numPC = 0;
end

%% Initialize for competitive phase
agconc_Epmask = epitopeMasking(agconc, abconc, Ka, param); 
                                        %Get effective concentration
numTcell = getNumTcells(param); % Get number of T cells
tspan = 0:param.dt:param.tmax;
tspan_summary = param.tspan_summary;
T = 0; % Total simulation time
activeGCnum = 1;
%% Competitive Phase Begins
for idx = 1:length(tspan) % At each time step

t = tspan(idx);
currentTcell = numTcell(idx,2);
activeGCnum = max(activeGCnum, ceil(currentTcell/1000));
borderTcell = numTcell(idx,1);
activeGCnum = min(param.M_GC, activeGCnum);
conc = [param.Ageff, 1]*agconc_Epmask;

%% Flux of naive B cells
% Get incoming naive B cells 
if t>0 
[incomingLogicIdx] = naiveFlux(naiveBcells, conc, param, borderTcell);
if any(any(incomingLogicIdx)) % If any naive B cell incoming
n_naive = size(naiveBcells,2);
entry_time_idx = [false(4*M_GC,n_naive);
                  incomingLogicIdx;
                  false(2*M_GC,n_naive)];
naiveBcells(entry_time_idx) = (t+0.5*param.dt); % time of entry
cnt = 0;
numIncoming = sum(incomingLogicIdx,2);
incoming = zeros(6,sum(numIncoming));

for k=1:M_GC
    incoming(:,cnt+1:cnt+numIncoming(k)) = ...
        naiveBcells(k+M_GC*[0:3,5:6], incomingLogicIdx(k,:));
    cnt = cnt+numIncoming(k);
end
GCidx = randi(activeGCnum, 1, sum(numIncoming));
% Add naive B cells to GCs
for i=1:sum(numIncoming)
    k = GCidx(i);
    if rand < param.outputprob
        if rand<param.outputpbfrac %&& any((0<t-param.dose_t) & (t-param.dose_t<5))
            plasmaBlasts(k:M_GC:(4*M_GC+k), (numPB(k)+1):(numPB(k)+32)) = ...
                repmat([incoming(1:4,i);t+0.5*param.dt],1,32);
            numPB(k) = numPB(k)+1;
        else
            memoryCells(k:M_GC:k+5*M_GC, numMC(k)+1) = ...
                [incoming(1:4,i);-1;t+0.5*param.dt];
            numMC(k) = numMC(k) + 1;
        end
    end
    copies = param.naiveprolnum; % Makes 4 copies
    % new GC B cells
    newIdx = (numGC(k)+1):(numGC(k)+copies);
    gcBcells(k:M_GC:(3*M_GC+k), newIdx) = ...
        repmat(incoming(1:4,i),1,copies);
    numGC(k) = numGC(k) + copies;
    % if the new GC B cells are from pre-existing memory cells
    if incoming(1,i)>2010 && any(incoming(5:6,i)>0)
       gcMutations(k,res2ind(newIdx,1:n_res,param))...
           = repmat(reshape(int2bit(incoming(5:6,i),n_res/2 ...
           ),n_res,1)',1,copies);
    end
    % Can check if the affinity calculation is correct, if wanted
   debugMutation = 0;
    if debugMutation
     for bid = newIdx
        checkMutationScheme(gcBcells, k, bid, naiveBcells,...
            mutations, gcMutations, M_GC, param);
     end
    end
end

end % end if
end
%% GC B cells birth and death
if any(any(numGC)) %If any GC B cell
[gcBirthIdx, gcDeathIdx] = birthDeath(gcBcells, conc, currentTcell, ...
    numGC, activeGCnum, param); %Get index of birth and death

% Birth
if any(any(gcBirthIdx))
    numBirths = sum(gcBirthIdx,2);
    dtCells = zeros(M_GC*5, max(numBirths));
        %Lineage, Target, Affinity-WT, Affinity-Var, Mutnum
    dtMuts = zeros(M_GC, max(numBirths));
        %Keep track of mutated residue, if a mutation occurs

    for k=1:M_GC %for each GC
        birthidx = find(gcBirthIdx(k,:));
        mutidx = zeros(size(birthidx)); %Keeps track of the actual 
        jj=1;                       % daughter cell indices in order to
                                    %copy over mutation state later

        if ~isempty(birthidx)
            for j=1:length(birthidx)
                    % Mutation can only happen to the daughter cell
                r1 = rand; r2 = rand;

                if r1<param.outputprob % positively selected cell 
                                        %differentiates into Mem or PC
                    if rand<param.outputpcfrac % becomes a plasma cell
                        plasmaCells(k:M_GC:(5*M_GC+k), ...
                                    numPC(k)+1) = ...
                            [gcBcells(k:M_GC:(4*M_GC+k), ...
                            birthidx(j));t+0.5*param.dt];
                        numPC(k) = numPC(k)+1;
                    else % becomes a memory cell
                        memoryCells(k:M_GC:k+5*M_GC, ...
                                    numMC(k)+1) = ...
                            [gcBcells(k:M_GC:k+4*M_GC, ...
                            birthidx(j));t+0.5*param.dt];
                        memoryCells(k+6*M_GC, numMC(k)+1) =...
                            bit2int(reshape(gcMutations(k,res2ind( ...
                            birthidx(j),1:n_res/2,param)),n_res/2,[]), ...
                            n_res/2);
                        memoryCells(k+7*M_GC, numMC(k)+1) =...
                            bit2int(reshape(gcMutations(k,res2ind( ...
                            birthidx(j),n_res/2+1:n_res,param)), ...
                            n_res/2,[]), n_res/2);
                        numMC(k) = numMC(k) + 1;
                    end
                    numBirths(k) = numBirths(k)-1;% no daughter cell
                    gcDeathIdx(k,birthidx(j)) = 1;% parent cell removed
                                           
                elseif r2<0.3 %Mutation leads to death
                    numBirths(k) = numBirths(k)-1; % no daughter cell
                elseif r2<0.8 || t<6 % No mutation until t=6
                    % Ref: - Jacob, Miller, Kelsoe 1992 Immunol. Cell. Bio
                    dtCells(k:M_GC:(4*M_GC+k), jj) = ...
                        gcBcells(k:M_GC:(4*M_GC+k), birthidx(j));
                    if t>6
                       dtCells(4*M_GC+k, jj) = ...
                           dtCells(4*M_GC+k, jj)+1; %Silent mutation
                    end
                    mutidx(jj) = birthidx(j);
                    jj = jj+1;

                else %Mutation that changes affinity
                    dtCells(k:M_GC:(4*M_GC+k), jj) = ...
                        gcBcells(k:M_GC:(4*M_GC+k), birthidx(j));
                    lineage = gcBcells(k, birthidx(j)); %get the lineage
                    mutres = randi(n_res); %Randomly pick residue
                    dE = mutations{1}(k, lineage, mutres); %WT affinity 
                    dEvar = mutations{2}(k, lineage, mutres); %Variant
                    already_mutated = gcMutations(k, ...
                                    res2ind(birthidx(j), mutres, param));
                    if already_mutated %Mutate from 1 to 0
                        dtMuts(k, jj) = -mutres;
                        dE = -dE;
                        dEvar = -dEvar;
                        dtCells(4*M_GC+k, jj) =  ...
                        dtCells(4*M_GC+k, jj)-100; %Affinity-changing
                    else
                        dtMuts(k, jj) = mutres;
                        dtCells(4*M_GC+k, jj) =  ...
                        dtCells(4*M_GC+k, jj)+100; %Affinity-changing
                    end

                    %Update the affinity + sanity check
                    dtCells(k+M_GC*2, jj) = ...
                        dtCells(k+M_GC*2, jj) - dE;
                    dtCells(k+M_GC*3, jj) = ...
                        dtCells(k+M_GC*3, jj) - dEvar;
                    if dE<-8 %no mutation can be this largely beneficial
                        save('debug.mat')
                        error('Error: beneficial mutation too large')
                    end
                    if any(dtCells(k+M_GC*[2;3], jj)>16)
                        save('debug.mat')
                        error('Error: affinity impossibly high')
                    end
                    mutidx(jj) = birthidx(j);
                    jj = jj+1;
                end %end if 
            end %end for

            % Keep track of the mutational state of each GC B cell
            if numBirths(k)>0
                newIdx = numGC(k)+1:numGC(k)+numBirths(k);
                gcBcells(k:M_GC:4*M_GC+k,newIdx)...
                    = dtCells(k:M_GC:4*M_GC+k,1:numBirths(k));
                mutDirection = (dtMuts(k,1:numBirths(k))>0)...
                    - (dtMuts(k,1:numBirths(k))<0); 
                            %1 if 0->1, -1 if 1->0, 0 if otherwise
                dtMuts(k,:) = abs(dtMuts(k,:)); 
                            %Index of mutated residue

                %Copy over the mutation state, then update it            
                for j=1:numBirths(k) 
                    gcMutations(k,res2ind(numGC(k)+j,1:n_res,param))...
                        = gcMutations(k,res2ind(mutidx(j),1:n_res,param));
                end
                mutated = find(mutDirection);
                resIdx = res2ind(newIdx(mutated),dtMuts(k,mutated),param);
                gcMutations(k, resIdx) = ...
                    gcMutations(k,resIdx)+mutDirection(mutated);
                debugMutation = 0;
                if debugMutation
                    for bid = newcellsidx
                       checkMutationScheme(gcBcells, k, bid, naiveBcells, ...
                           mutations, gcBcellsMutation, M_GC, param);
                    end
                end
                numGC(k) = numGC(k) + numBirths(k);
            end
        end
    end
end

%Deaths
for k=1:M_GC %Death or exit as memory
    deathidx = find(gcDeathIdx(k,:));
    gcBcells(k:M_GC:k+4*M_GC, deathidx) = 0;
    for j = 1:length(deathidx)
        gcMutations(k, res2ind(deathidx(j),1:n_res,param))=0;
    end

    %Cleanup of the arrays to remove the 0s in betwen and align to left
    if numGC(k)>N_GC_MAX*0.95
        liveBcells = find(gcBcells(k,:));
        gcBcells(k:M_GC:k+4*M_GC, 1:length(liveBcells)) = ...
            gcBcells(k:M_GC:k+4*M_GC, liveBcells);
        gcBcells(k:M_GC:k+4*M_GC, length(liveBcells)+1:end) = 0;
        numGC(k) = length(liveBcells);
        gcMutations(k, res2ind(1:numGC(k), 1:n_res, param)) = ...
            gcMutations(k, res2ind(liveBcells, 1:n_res, param));
        gcMutations(k, res2ind(numGC(k)+1, 1, param):end)=0;
        if numGC(k)>N_GC_MAX*0.95 %If array is almost filled, then expand
            gcBcells = [gcBcells,flatten2D(zeros( ...
                M_GC, param.N_GC_MAX/2, param.gcfieldnum))];
            gcMutations=[gcMutations, zeros(M_GC, param.N_GC_MAX/2*n_res)];
            N_GC_MAX = N_GC_MAX+param.N_GC_MAX/2;
        end
    end
end
end

%% Output cells birth and death
% EGC Birth
if any(((t-param.dose_t(2:end))<6) & ((t-param.dose_t(2:end))>-0.005)) ||...
    (param.pSER>0 && (t-param.dose_t(end))>-0.005) && (t-param.dose_t(end)<param.pSER_ag_release_time+6)
    if any(abs(t-param.dose_t(2:end))<0.005)
        buildup3D = @(A,d3) permute(reshape(A',size(A',1),[],d3),[2,1,3]);
        existingMemory = buildup3D(memoryCells, param.memfieldnum);
        flatten1D = @(A) reshape(permute(A,[3,2,1]), size(A,3),[]);
        existingMemory = alignLeft(flatten1D(existingMemory), 1);
        memoryCellsEGC(:,1:size(existingMemory,2)) = existingMemory;
        numMCEGC = numMCEGC + size(existingMemory,2);
        memoryCells(:) = 0;
        numMC(:) = 0;
    end

    [memBirthIdx, plasmaIdx] = ...
    birthDeathEGC(memoryCellsEGC(:,1:numMCEGC), conc, param);
    if any(plasmaIdx)
    newIdx = numPCEGC+1:numPCEGC+length(plasmaIdx);
    plasmaCellsEGC(1:5,newIdx) = memoryCellsEGC(1:5,plasmaIdx);
    plasmaCellsEGC(6, newIdx) = t;
    numPCEGC = numPCEGC+length(plasmaIdx);
    end
    if any(memBirthIdx)
    newIdx = numMCEGC+1:numMCEGC+length(memBirthIdx);
    memoryCellsEGC(:,newIdx) = memoryCellsEGC(:,memBirthIdx);
    numMCEGC = numMCEGC + length(memBirthIdx);
    end
end

% EGC Death
[pcDeathIdx,egcPCDeathIdx] =...
  memPCDeath(plasmaCells,plasmaCellsEGC,param);
dead.plasmaCellsEGC(:,dead.numPC+1:dead.numPC+sum(egcPCDeathIdx)) = ...
    plasmaCellsEGC(:,egcPCDeathIdx);
dead.numPC = dead.numPC + sum(egcPCDeathIdx);
temp = plasmaCells';
deadPCs = reshape(temp(repmat(pcDeathIdx,param.pcfieldnum,1)'), ...
                  [],param.pcfieldnum)';
dead.plasmaCells(:,dead.PCnum+1:dead.PCnum+size(deadPCs,2)) = deadPCs;
dead.PCnum = dead.PCnum + size(deadPCs,2);
if dead.numPC>0.95*length(dead.plasmaCellsEGC)
    dead.plasmaCellsEGC = [dead.plasmaCellsEGC, ...
                          zeros(size(dead.plasmaCellsEGC))];
end
if dead.PCnum>0.95*length(dead.plasmaCells)
    dead.plasmaCells = [dead.plasmaCells, zeros(size(dead.plasmaCells))];
end
plasmaCellsEGC(:,egcPCDeathIdx) = 0;
plasmaCells(repmat(pcDeathIdx,param.pcfieldnum,1)) = 0;

%Array scaling
if any(numMC>(0.95*size(memoryCells,2)))
    memoryCells = [memoryCells, zeros(size(memoryCells))];
end
if numMCEGC > 0.95*size(memoryCellsEGC,2)
    memoryCellsEGC = [memoryCellsEGC, zeros(size(memoryCellsEGC))];
end
if numPCEGC > 0.95*size(plasmaCellsEGC,2)
    plasmaCellsEGC = alignLeft(plasmaCellsEGC, 0);
    numPCEGC = find(plasmaCellsEGC(1,:),1,'last');
    if numPCEGC > 0.95*size(plasmaCellsEGC,2)
        plasmaCellsEGC = [plasmaCellsEGC, zeros(size(plasmaCellsEGC))];
    end
end

T = T + toc(tstart);
%% Update the concentrations
if param.numshot == 100
    %This triggers the condition for manual Ag concentration control
    param.k_deposit = 0; %IC deposition is suppressed
    param.dose_t = [];
    param.dose = [];
    param.d_IC = 0;
end

[agconc, abconc, Ka, Ka_var] = ...
  updateConcentrations(agconc, abconc, Ka, Ka_var, ...
  plasmaBlasts(:,1:max(max(numPB),1)), plasmaCells(:,1:max(max(numPC),1)), ...
  plasmaCellsEGC(:,1:max(numPCEGC,1)), t, param);
agconc_Epmask = epitopeMasking(agconc, abconc, Ka, param);
if param.numshot == 100
    agconc_Epmask = setConcManual(t, param.k, param.T, param.C0, param.masking);
    agconc(2) = agconc_Epmask(2,1);
end


%% Store results intermittently
if ismember(t, tspan_summary)
    storeidx = find(t==tspan_summary);
    result = updateHistory(result, gcBcells, plasmaCells, memoryCells,...
        plasmaCellsEGC, memoryCellsEGC, agconc, abconc, Ka, Ka_var,...
        agconc_Epmask, tspan_summary, storeidx, param);
end
if ismember(t, 1:1:180)
   fprintf(['t=',num2str(t),'\n']);
   toc(tstart)
end

end % END of competitive phase
T

%% Cleanup and store result
for k=1:M_GC
liveBcells = find(gcBcells(k,:));
gcBcells(k:M_GC:k+4*M_GC, 1:length(liveBcells)) = ...
    gcBcells(k:M_GC:k+4*M_GC, liveBcells);
gcBcells(k:M_GC:k+4*M_GC, length(liveBcells)+1:end) = 0;
numGC(k) = length(liveBcells);
gcMutations(k, res2ind(1:numGC(k), 1:n_res, param)) = ...
    gcMutations(k, res2ind(liveBcells, 1:n_res, param));
gcMutations(k, res2ind(numGC(k)+1, 1, param):end)=0;
plasmaBlasts(k:M_GC:k+(param.pcfieldnum-2)*M_GC, :) = ...
    alignLeft(plasmaBlasts(k:M_GC:k+(param.pcfieldnum-2)*M_GC, :),0);
plasmaCells(k:M_GC:k+(param.pcfieldnum-1)*M_GC, :) = ...
    alignLeft(plasmaCells(k:M_GC:k+(param.pcfieldnum-1)*M_GC, :),0);
memoryCells(k:M_GC:k+(param.memfieldnum-1)*M_GC, :) = ...
    alignLeft(memoryCells(k:M_GC:k+(param.memfieldnum-1)*M_GC, :),0);
end

plasmaCellsEGC = alignLeft(plasmaCellsEGC, 0);
memoryCellsEGC = alignLeft(memoryCellsEGC, 0);

%Store final results
buildup3D = @(A,d3) permute(reshape(A',size(A',1),[],d3),[2,1,3]);
result.output.finalpb = buildup3D(plasmaBlasts, param.pcfieldnum-1);
result.output.finalpc = buildup3D(plasmaCells, param.pcfieldnum);
result.output.finalmem = buildup3D(memoryCells, param.memfieldnum);
result.naive = buildup3D(naiveBcells, param.naivefieldnum);
result.gc.finalgc = buildup3D(gcBcells, param.gcfieldnum);
result.memoryCellsEGC = memoryCellsEGC;
result.plasmaCellsEGC = plasmaCellsEGC;
result.memoryCellsEGC = memoryCellsEGC;
result.dead = dead;

idx = 1+4*(0:7:param.tmax);
result.gc.numbylineage_top = zeros(param.M_GC,length(param.tspan_summary),10);
for i=1:param.M_GC
    topidx = find(result.naive(i,:,3)==max(result.naive(i,:,3)),1,'last');
    topidx = topidx-9:topidx;
    result.gc.numbylineage_top(i,:,:) = ...
        result.gc.numbylineage(i,:,topidx);
end
result.gc.numbylineage = result.gc.numbylineage(:,idx,:);
[result.output.finalmem, result.output.memMutations] = ...
    getMemoryMutations(result.output.finalmem, mutations, param);

%Reshape back to 3D arrays
if saveresult
    dirnm = fileparts(fnm{param.vaxnum});
    if ~isfolder(dirnm)
        mkdir(dirnm)
    end
    save(fnm{param.vaxnum}, 'result')
end
toc(tstart)

end

%#################################################################
%% subfunctions
function [memory, mutationsCompact] = ...
    getMemoryMutations(memory, mutations, param)
%-----------------------------------------------------------------
% Summary:
%  Returns the mutation sizes of the memory cells. This enables only the
%  mutation sizes of relevant memory cells to be saved, for storage space.
%  Then, this information can be accessed upon next vaccination.

% Outputs:
%  memory: Array of GC-derived memory cells. 
%  mutationsCompact:2 x (max number of memory cells in GC) x 80 array

% Inputs:
%  memory: Array of GC-derived memory cells. 
%  mutations: 1x2 cell array of mutations
%  param: parameter struct
%-----------------------------------------------------------------
[m1,m2,m3] = size(memory);
uniqueCloneNum = 0;
mutationsCompact = zeros(2, m2, param.n_res);
for k=1:m1
   uniqueClones = setdiff(unique(squeeze(memory(k,:,1))),[0]);
   for j=1:length(uniqueClones)
      uniqueCloneNum = uniqueCloneNum + 1;
      memory(k,squeeze(memory(k,:,1)==uniqueClones(j)),9) = uniqueCloneNum;
      for i=1:2
          mutationsCompact(i,uniqueCloneNum,:) = ....
              mutations{i}(k,uniqueClones(j),:);
      end
   end
end
mutationsCompact = mutationsCompact(:,1:uniqueCloneNum,:);
end



function [pcDeathIdx, egcPCDeathIdx] =...
    memPCDeath(plasmaCells, plasmaCellsEGC, param)
%-----------------------------------------------------------------
% Summary: Find indices of PCs that will undergo apoptosis
% Inputs: GC and EGC-derived PC arrays, parameter struct
% Outputs: 
%  pcDeathIdx: M_GC x (GC-derived PC max capacity) logical indices
%  egcPCDeathIdx: 1 x (EGC-derived PC max capacity) logical indices
%-----------------------------------------------------------------
M_GC = param.M_GC;
pcDeathIdx = plasmaCells(1:M_GC,:)>0 & ...
    rand(size(plasmaCells(1:M_GC,:)))<param.d_pc*param.dt;
egcPCDeathIdx = plasmaCellsEGC(1,:)>0 & ...
    rand(size(plasmaCellsEGC(1,:)))<param.d_pc*param.dt;
end

function [memBirthIdx, plasmaIdx] = ...
         birthDeathEGC(memoryCellsEGC, conc, param)
%-----------------------------------------------------------------
% Summary:
%  Finds the indices of EGC memory B cells that will give birth to either
%  new memory cells or to plasma cells

% Outputs:
%  memBirthIdx, plasmaIdx: Indices (not logical). Length sums up to the
%  number of EGC memory cells. 

% Inputs: EGC-derived memory cells, concentration, parameter struct
%-----------------------------------------------------------------
concarray = zeros(size(memoryCellsEGC(2,:)));
for ep = 1:param.n_ep
   concarray = concarray + conc(ep)*(memoryCellsEGC(2,:)==ep);
end
activation_signal = ((concarray/param.C0).^0.5.*...
    (10.^(min(memoryCellsEGC(3,:),10)-param.f0))).^param.w2;
if param.w1 > 0 
    activation_signal = ((param.w1+1)*activation_signal)./ ...
                        (param.w1+activation_signal);
end
activation = min(activation_signal,1) > rand(size(activation_signal));
% activation = ones(size(activation_signal));
n_activated = sum(activation,2);
if any(n_activated) %at least one B cell is intrinsically activated
    %% Extrinsic activation by T cells
    % define competitiveness of each precursor
    activated_fitness = activation.*(activation_signal);
    avgfitness = mean(activated_fitness(activated_fitness>0));
    helpamount = (param.M_GC*param.numTmax./n_activated).* ...
                  activated_fitness./avgfitness;
    beta = param.betamaxEGC*(helpamount./(1+helpamount));
    beta(isnan(beta))=0;
    temp = find(beta);
    temp = temp(rand(size(temp))<beta(temp)*param.dt);
    r = rand(size(temp))<0.6;
    plasmaIdx = temp(r);
    memBirthIdx = temp(~r);
else
    memBirthIdx = [];
    plasmaIdx = [];
end
end

function [gcBirthIdx, gcDeathIdx] =...
    birthDeath(gcBcells, conc, currentTcell, currentGCnum, activeGCnum, param)
%-----------------------------------------------------------------
% Summary:
%  Finds the indices of GC B cells that are positively selected and that
%  will undergo apoptosis

% Outputs:
%  gcBirthIdx, gcDeathIdx: M_GC x length of GC arrays. Logical indices

% Inputs:
%  GC B cells, concentration, T cell number, current number of B cells in
%  GC, parameter struct
%-----------------------------------------------------------------
M_GC = param.M_GC;
concarray = zeros(size(gcBcells(M_GC+1:2*M_GC,:)));
for ep = 1:param.n_ep
   concarray = concarray + conc(ep)*(gcBcells(M_GC+1:2*M_GC,:)==ep);
end
activation_signal = ((concarray/param.C0).^0.5.*...
    (10.^(min(gcBcells(M_GC*2+1:M_GC*3,:),10)-param.f0))).^param.w2;
if param.w1 > 0 
    activation_signal = ((param.w1+1)*activation_signal)./ ...
                        (param.w1+activation_signal);
end
activation = min(activation_signal,1) > rand(size(activation_signal));
activation(activeGCnum+1:end, :) = 0;
n_activated = sum(activation,2);

if any(n_activated) %at least one B cell is intrinsically activated
    %% Extrinsic activation by T cells
    % define competitiveness of each precursor
    activated_fitness = activation.*(activation_signal);
    avgfitness = zeros(M_GC,1);
    for k=1:M_GC
        avgfitness(k) = mean(activated_fitness(k,...
                         activated_fitness(k,:)>0));
    end
    helpamount = ((currentTcell/activeGCnum)./n_activated).* ...
                    activated_fitness./avgfitness;
    beta = param.betamax*(helpamount./(1+helpamount));
    beta(isnan(beta))=0;
    temp = find(beta);
    temp = temp(rand(size(temp))<beta(temp)*param.dt);
    gcBirthIdx = false(size(beta));
    gcBirthIdx(temp) = 1;
else
    gcBirthIdx = [];
end
gcDeathIdx = false(size(activation));
for k=1:M_GC
    live = gcBcells(k, 1:currentGCnum(k))~=0;
%     mu = param.mu*(sum(live)/param.N_GC_MAX);
    mu = param.mu;
    gcDeathIdx(k, 1:currentGCnum(k)) = rand(1,currentGCnum(k)) ...
                                        <mu*param.dt & live;
end
end

function [incomingnaive] = naiveFlux(naiveBcells, conc, param, borderTcell)
%-----------------------------------------------------------------
% Summary:
%  Finds the logical indices of naive B cells that will enter GC

% Outputs:
%  incomingnaive: M_GC x length of naive B cells array of logical indices

% Inputs:
%  Array of naive B cells, concentration, parameter struct
%-----------------------------------------------------------------
M_GC = param.M_GC;
concarray = zeros(size(naiveBcells(M_GC+1:2*M_GC,:)));
for ep = 1:param.n_ep
   concarray = concarray + conc(ep)*(naiveBcells(M_GC+1:2*M_GC,:) ...
       ==ep).*(naiveBcells(M_GC*4+1:5*M_GC,:)==0);
end
activation_signal = ((concarray/param.C0).^0.5.*...
    (10.^(naiveBcells(M_GC*2+1:M_GC*3,:)-param.f0))).^param.w2;
% activation_signal(activation_signal<0.2) = 0;
if param.w1 > 0 
    activation_signal = ((param.w1+1)*activation_signal)./ ...
    (param.w1+activation_signal);
end
% encounter = 0.01 > rand(size(activation_signal));
activation = min(activation_signal,1) > rand(size(activation_signal));
% activation = activation & encounter;
n_activated = sum(sum(activation));

if any(n_activated) %at least one B cell is intrinsically activated
    %% Extrinsic activation by T cells
    % define competitiveness of each precursor
    activated_fitness = activation.*(activation_signal);
    avgfitness = zeros(M_GC,1);
    for k=1:M_GC
        avgfitness(k) = mean(activated_fitness( ...
                            k, activated_fitness(k,:)>0));
    end
    helpamount = (borderTcell./n_activated).*activated_fitness./avgfitness;
    lambda = param.lambdamax*(helpamount./(1+helpamount));
    incomingnaive = activation & rand(size(activation))<lambda*param.dt;
    
    %% Differentiation
else
    incomingnaive = [];
end
end



function [memoryToGC, mutations, memoryToEGC] =...
    splitMemory(memoryCellsEGC, memoryMutations, mutations, param)
%-----------------------------------------------------------------
% Summary:
%   Takes in the pre-existing memory cells and their fitness landscapes,
%   then selects a pre-define fraction of the memory B cells to be added to
%   naive B cells for potential re-activation

% Outputs:
%   memoryToGC: (M_GC*9) x MemoryReentryMax array.
%               Memory cells that will be added to naive B cells
%   mutations: 1 x 2 cell array of mutations
%   memoryToEGC: 9 x N, N is the number of memory cells that goes to EGC

% Inputs:
%   memoryCellsEGC: EGC-derived memory cells 
%   memoryMutations: 2 x N x 80 array. Mutation sizes of all unique memory
%                    cells. Dim1 - target; Dim2 - Number of unique memory
%                    Dim3 - Residues
%   mutations: 1x2 cell of mutation sizes
%   param: parameter struct
%-----------------------------------------------------------------
%% Function main 
rng(1)

% Select the memory cells to be added to re-activation pool
N_mem = find(memoryCellsEGC(1,:),1,'last');
toGC = binornd(ones(1,N_mem), param.memToGCFrac);
memoryToGC = memoryCellsEGC(:,toGC==1);
memoryToGC = alignLeft(memoryToGC, 1);
memoryCellsEGC(:,toGC==1)=0;
memoryToEGC = alignLeft(memoryCellsEGC, 0);

% Distribute the memory cells over the secondary GCs
N_memToGC = size(memoryToGC,2);
memPoolGC = zeros(param.M_GC, param.MemoryReentryMax, param.naivefieldnum);
memNumGC = zeros(param.M_GC,1);
for i=1:N_memToGC % Iterative over selected memory cells
   k = randi(param.M_GC); % Randomly select GC and assign it there
   memNumGC(k) = memNumGC(k)+1;
   memPoolGC(k, memNumGC(k), 1) = memNumGC(k)+param.NaiveMax; % lineage
   memPoolGC(k, memNumGC(k), 2:4) = memoryToGC(2:4, i); %target, affinities
   memPoolGC(k, memNumGC(k), 6:7) = memoryToGC(7:8, i); %mutation states
   for j=1:2 % Copy over mutation states
       mutations{j}(k, memNumGC(k)+param.NaiveMax, :) =...
           memoryMutations(j, memoryToGC(9,i), :);
   end
   
end
flatten2D = @(A) reshape(permute(A,[2,1,3]), size(A,2),[])';
memoryToGC = flatten2D(memPoolGC);

end
