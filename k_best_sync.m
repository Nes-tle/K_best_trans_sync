function [Results] = k_best_sync(pairmatches, Para)
% Perform k best synchronization on pair-wise matches between pairs of
% range scans
% 'pairmatches' store the pairwise matches:
%  pairmatches{id}.sId : the source scan index of this pairwise match 
%  pairmatches{id}.tId : the target scan index of this pairwise match
%  pairmatches{id}.R_opt: the rotational component of this pairwise match
%  pairmatches{id}.t_opt: the translation component of this pairwise match
%  Para.numIters : number of iterations of the alternating procedures
%  Para.maxK: the maximum size of the K, the value of K is determined
%  automatically
% Result: Result is a cell of dimension numScans
%  
[Neighbors, rootId, numScans] = compute_neighbor_matches(pairmatches);
% Initialize the poses
[K, grid_size] = shallow_analysis(Neighbors, rootId, numScans, Para);
Para.K = K;
Para.grid_size = grid_size;
% 
Results = cell(1, numScans);
for sId = 1 : numScans
    Results{sId}.R_opts = [];
    Results{sId}.t_opts = [];
end
Results{rootId}.R_opts{1} = eye(3);
Results{rootId}.t_opts{1} = zeros(3,1);
Results{rootId}.weights = 1;
for iteration = 1 : Para.maxNumIters
    Results_temp = propagate_poses(Results, Neighbors);
    for scanId = 1 : numScans
        Results_temp{scanId} = clustering(...
            Results_temp{scanId}, Para.grid_size, Para.K);
    end
    if check_convergence(Results, Results_temp, rootId, Para) == 1
        Results = Results_temp;
        break;
    end
    Results = Results_temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We perform a few iterations to analyze the data-dependent parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, grid_size] = shallow_analysis(Neighbors, rootId, numScans, Para)
%
Results = cell(1, numScans);
for sId = 1 : numScans
    Results{sId}.R_opts = [];
    Results{sId}.t_opts = [];
end
Results{rootId}.R_opts{1} = eye(3);
Results{rootId}.t_opts{1} = zeros(3,1);
Results{rootId}.weights = 1;
% Analyze the propagted poses
Results = propagate_poses(Results, Neighbors);
Results = propagate_poses(Results, Neighbors);
densities = zeros(1, numScans);
for sId =  1 : numScans
    densities(sId) = extract_density(Results{sId});
end
grid_size = median(densities);
votsforK = zeros(1, numScans);
for sId = 1 : numScans
    Result_clus = clustering(Results{sId}, Para.grid_size_rot, Para.maxK);
    ws = Result_clus.weights;
    curve = ws(1:(length(ws)-2)) - 2*ws(2:(length(ws)-1)) + ws(3:length(ws));
    [~,off] = max(curve);
    votsforK(sId) = off;
end
% calculate the majority vote
J = sparse(1:numScans, votsforK, ones(1,numScans));
[~, off] = max(sum(J));
K = off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decide whether two groups of solutions are consistent or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flag] = check_convergence(Results, Results_temp, root_Id, Para)
% Determine if we want to stop the recursion or not
% Check the condition one
% When the number of solutions per scan is not fully reached, then the recursion 
% has not converged of course. 
flag = 1;
for sId = 1 : length(Results)
    if min(length(Results{sId}.R_opts),length(Results_temp{sId}.R_opts)) < Para.K
        flag = 0;
        break;
    end
end
if flag == 0
    return;
end
numScans = length(Results);

% Otherwise, compare the resulting transformations 
Results_aug = cell(1, numScans);
Results_temp_aug = cell(1, numScans);
for sId = 1 : numScans
    Results_aug{sId} = zeros(12, Para.K*Para.K);
    Results_temp_aug{sId} = zeros(12, Para.K*Para.K);
end
for sId = 1 : numScans
    for i = 1 : Para.K
        Rs = Results{sId}.R_opts{i};
        ts = Results{sId}.t_opts{i};
        Rs_temp = Results_temp{sId}.R_opts{i};
        ts_temp = Results_temp{sId}.t_opts{i};
        for j = 1 : Para.K
            Rs_root = Results{root_Id}.R_opts{j};
            ts_root = Results{root_Id}.t_opts{j};
            Rs_root_temp = Results_temp{root_Id}.R_opts{j};
            ts_root_temp = Results_temp{root_Id}.t_opts{j};
            off = (i-1)*Para.K + j;
            R = Rs*Rs_root';
            t = ts - R*ts_root;
            Results_aug{sId}(1:9, off) = reshape(R, [9,1]);
            Results_aug{sId}(10:12,off) = t;
            % Normalize based on the root scan
            R_temp = Rs_temp*Rs_root_temp';
            t_temp = ts_temp - R_temp*ts_root_temp;
            Results_temp_aug{sId}(1:9, off) = reshape(R_temp,[9,1]);
            Results_temp_aug{sId}(10:12,off) = t_temp;
        end
    end
end
mean_error = 0;
for sId =  1 : numScans
    error = distance(Results_aug{sId}, Results_temp_aug{sId});
    mean_error = mean_error + error;
end
mean_error = mean_error/numScans;
fprintf('%f\n', mean_error);
if mean_error < Para.eps
    flag = 1;
else
    flag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the sampling density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [density] = extract_density(Result)
%
num = length(Result.R_opts);
data = zeros(12, num);
for id = 1 : num
    data(1:9, id) = reshape(Result.R_opts{id}, [9,1]);
    data(10:12,id) = Result.t_opts{id};
end
[IDX, DIS] = knnsearch(data', data', 'k', min(10, num));
density = min(median(DIS'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum distance between two point clouds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [error] = distance(points1, points2)
%
num1 = size(points1, 2);
num2 = size(points2, 2);
ids1 = kron(ones(1,num2), 1:num1);
ids2 = kron(1:num2, ones(1,num1));
dif = points1(:, ids1) - points2(:, ids2);
dif = ones(1, size(points1,1))*(dif.*dif);
dif = sqrt(reshape(dif, [num1, num2]));
error = (mean(min(dif)) + mean(min(dif')))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform clustering for each pose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Result_clus] = clustering(Result, grid_size, maxnumClus)
if length(Result.R_opts) <= maxnumClus
    Result_clus = Result;
    return;
end
data = zeros(12, length(Result.R_opts));
weights = zeros(1, length(Result.R_opts));
for id = 1 : length(Result.R_opts)
    data(1:9, id) = reshape(Result.R_opts{id}, [9,1]);
    data(10:12, id) = Result.t_opts{id};
    weights(id) = Result.weights(id);
end
% This computation can be accelerated
num = size(data,2);
ids1 = kron(ones(1,num),1:num);
ids2 = kron(1:num, ones(1,num));
sigma = grid_size;
for outer = 1 : 20
    disMat = data(:,ids1) - data(:,ids2);
    disMat = sum(disMat.*disMat);
    weightMat = exp(-disMat/2/sigma/sigma);
    weightMat = reshape(weightMat, [num,num]);
    weightMat = weightMat.*(weights'*ones(1,num));
    data = (data*weightMat)./(ones(12,1)*sum(weightMat));
end
%
clusterIds = zeros(1, size(data,2));
numClus = 0;
for id = 1 : size(data, 2)
    if clusterIds(id) == 0
        numClus = numClus + 1;
        clusterIds(id) = numClus;
        for id2 = (id+1) : size(data,2)
            dif = data(:,id) - data(:, id2);
            if norm(dif) < 1e-2
                clusterIds(id2) = numClus;
            end
        end
    end
end
J = sparse(1:length(clusterIds), clusterIds, ones(1,length(clusterIds)),...
 length(clusterIds), numClus);
clusters = zeros(13,numClus);
for id = 1 : numClus
    ids = find(J(:,id));
    v = ones(1,length(ids))*(data(:,ids)'.*(weights(ids)'*ones(1,12)))/sum(weights(ids));
    clusters(1:12,id) = v;
    clusters(13,id) = sum(weights(ids));
end
if numClus > maxnumClus
    [s, order] = sort(-clusters(13,:));
    clusters = clusters(:, order(1:numClus));
    numClus = maxnumClus;
end
clusters(13,:) = clusters(13,:)/sum(clusters(13,:));
for id = 1 : numClus
    v = clusters(1:12, id);
    R_opt = reshape(v(1:9), [3,3]);
    t_opt = v(10:12);
    [u,v,w] = svd(R_opt);
    R_opt = u*w';
    Result_clus.R_opts{id} = R_opt;
    Result_clus.t_opts{id} = t_opt;
    Result_clus.weights(id) = clusters(13, id);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform one step of propagation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Results_out] = propagate_poses(Results_in, Neighbors)
%
numScans = length(Results_in);
Results_out = cell(1, numScans);
for sId = 1 : numScans
    Results_out{sId}.R_opts = [];
    Results_out{sId}.t_opts = [];
    Results_out{sId}.weights = [];
end
%
for sId = 1 : numScans
    Result_in = Results_in{sId};
    for k = 1 : length(Result_in.R_opts)
        R_s = Result_in.R_opts{k};
        t_s = Result_in.t_opts{k};
        w = Result_in.weights(k);
        for j = 1 : length(Neighbors{sId}.tIds)
            tId = Neighbors{sId}.tIds(j);
            R_st = Neighbors{sId}.R_opts{j};
            t_st = Neighbors{sId}.t_opts{j};
            R_t = R_st*R_s;
            t_t = R_st*t_s + t_st;
            off = length(Results_out{tId}.R_opts);
            Results_out{tId}.R_opts{off+1} = R_t;
            Results_out{tId}.t_opts{off+1} = t_t;
            Results_out{tId}.weights(off+1) = w;
        end
    end
end

function [Neighbors, rootId, numScans] = compute_neighbor_matches(pairmatches)
% Compute the adjacency matrix
numEdges = size(pairmatches, 2);
edges = zeros(2, numEdges);
for eId = 1 : numEdges
    edges(1, eId) = pairmatches{eId}.sId;
    edges(2, eId) = pairmatches{eId}.tId;
end
numScans = max(max(edges));
A = sparse(edges(1,:), edges(2,:), ones(1,numEdges), numScans, numScans);
A = max(A, A');
vertexDegrees = full(sum(A));
% root is set as the scan with the maximum degree
[s, rootId] = max(vertexDegrees);
%
Neighbors = cell(1, numScans);
for scanId = 1 : numScans
    Neighbors{scanId}.tIds = [];
    Neighbors{scanId}.R_opts = [];
    Neighbors{scanId}.t_opts = [];
end
for pairId = 1 : length(pairmatches)
    pair = pairmatches{pairId};
    % Push from pair.sId to pair.tId
    Neighbors{pair.sId}.tIds = [Neighbors{pair.sId}.tIds, pair.tId];
    off = length(Neighbors{pair.sId}.tIds);
    Neighbors{pair.sId}.R_opts{off} = pair.R_opt;
    Neighbors{pair.sId}.t_opts{off} = pair.t_opt;
    % Push from pair.tId to pair.sId
    Neighbors{pair.tId}.tIds = [Neighbors{pair.tId}.tIds, pair.sId];
    off = length(Neighbors{pair.tId}.tIds);
    Neighbors{pair.tId}.R_opts{off} = pair.R_opt';
    Neighbors{pair.tId}.t_opts{off} = -pair.R_opt'*pair.t_opt;
end
% 
%     sub_clusters = rotation_based_clustering(...
%         Result.R_opts(clusters{id}), Result.weights(clusters{id}), Para);
%     for i = 1 : length(sub_clusters)
%         ids = clusters{id};
%         ids = ids(sub_clusters{i});
%         R_opt = zeros(3,3);
%         t_opt = zeros(3,1);
%         for j = 1 : length(ids)
%             R_opt = R_opt + Result.R_opts{ids(j)};
%             t_opt = t_opt + Result.t_opts{ids(j)};
%         end
%         R_opt = R_opt/length(ids);
%         t_opt = t_opt/length(ids);
%         [u,v,w] = svd(R_opt);
%         R_opt = u*w';
%         numClus = numClus + 1;
%         Result_clus.R_opts{numClus} = R_opt;
%         Result_clus.t_opts{numClus} = t_opt;
%         Result_clus.weights(numClus) = length(ids);
%     end
% end
% [s,order] = sort(-Result_clus.weights);
% order = order(1:Para.K);
% Result_clus.R_opts = Result_clus.R_opts(order);
% Result_clus.t_opts = Result_clus.t_opts(order);
% Result_clus.weights = Result_clus.weights(order);
% h = 10;
% % % 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Rotation-based clustering
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % function [clusters] = rotation_based_clustering(R_opts, weights, Para)
% % % if length(R_opts) == 0
% % %     clusters = [];
% %     return;
% % end
% % if length(R_opts) == 1
% %     clusters{1} = [1];
% %     return;
% % end
% % num = length(R_opts);
% % data = zeros(9, num);
% % for id = 1 : num
% %     data(:,id) = reshape(R_opts{id}, [9,1]);
% % end
% % ids1 = kron(ones(1,num),1:num);
% % ids2 = kron(1:num, ones(1,num));
% % sigma = Para.grid_size_rot;
% % for outer = 1 : 10
% %     disMat = data(:,ids1) - data(:,ids2);
% %     disMat = sum(disMat.*disMat);
% %     weightMat = exp(-disMat/2/sigma/sigma);
% %     weightMat = reshape(weightMat, [num,num]);
% %     weightMat = weightMat.*(weights'*ones(1,num));
% %     data = (data*weightMat)./(ones(9,1)*sum(weightMat));
% % end
% % clusterIds = zeros(1, size(data,2));
% % numClus = 0;
% % for id = 1 : size(data, 2)
% %     if clusterIds(id) == 0
% %         numClus = numClus + 1;
% %         clusterIds(id) = numClus;
% %         for id2 = (id+1) : size(data,2)
% %             dif = data(:,id) - data(:, id2);
% %             if norm(dif) < 1e-2
% %                 clusterIds(id2) = numClus;
% %             end
% %         end
% %     end
% % end
% % J = sparse(1:length(clusterIds), clusterIds, ones(1,length(clusterIds)),...
% %  length(clusterIds), numClus);
% % clusters = cell(1,numClus);
% % for id = 1 : numClus
% %     clusters{id} = find(J(:,id))';
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Translation-based clustering
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [clusters] = translation_based_clustering(t_opts, Para)
% % %
% % t_min = min(t_opts')';
% % t_max = max(t_opts')';
% % center = (t_max+t_min)/2;
% % dimx = floor((t_max(1)-t_min(1))/Para.grid_size_trans) + 3;
% % dimy = floor((t_max(2)-t_min(2))/Para.grid_size_trans) + 3;
% % dimz = floor((t_max(3)-t_min(3))/Para.grid_size_trans) + 3;
% % lower = [center(1) - dimx*Para.grid_size_trans/2;
% %     center(2) - dimy*Para.grid_size_trans/2;
% %     center(3) - dimz*Para.grid_size_trans/2];
% % upper = [center(1) + dimx*Para.grid_size_trans/2;
% %     center(2) + dimy*Para.grid_size_trans/2;
% %     center(3) + dimz*Para.grid_size_trans/2];
% % %
% % idsx = floor((t_opts(1,:) - lower(1))/Para.grid_size_trans);
% % idsy = floor((t_opts(2,:) - lower(2))/Para.grid_size_trans);
% % idsz = floor((t_opts(3,:) - lower(3))/Para.grid_size_trans);
% % num = size(t_opts, 2);
% % RowsJ = ones(27,1)*(1:num);
% % ColsJ = zeros(27, num);
% % ValsJ = ones(27, num);
% % for i = -1 : 1
% %     for j = -1 : 1
% %         for k = -1 : 1
% %             off = (i+1)*9 + (j+1)*3 + k + 2;
% %             ColsJ(off,:) = (idsx+i)*dimy*dimz + (idsy+j)*dimz + idsz + k + 1;
% %         end
% %     end
% % end
% % J = sparse(RowsJ, ColsJ, ValsJ, num, dimx*dimy*dimz);
% % colIds = find(sum(J));
% % J = J(:, colIds);
% % numClus = length(colIds);
% % for id = 1 : numClus
% %     clusters{id} = find(J(:,id));
% % end
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Perform clustering using a generic approach
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [Result_clus] = distance_clustering2(Result, Para)%
% % if length(Result.R_opts) <= Para.K
% %     Result_clus = Result;
% %     return;
% % end
% % num = length(Result.R_opts);
% % data_rot = zeros(9, num);
% % data_trans = zeros(3, num);
% % for id = 1 : num
% %     data_rot(:,id) = reshape(Result.R_opts{id}, [9,1]);
% %     data_trans(:,id) = Result.t_opts{id};
% % end
% % ids1 = kron(1:num, ones(1,num));
% % ids2 = kron(ones(1,num),1:num);
% % dist_rot = data_rot(:,ids1) -data_rot(:,ids2);
% % dist_rot = sqrt(ones(1,9)*(dist_rot.*dist_rot));
% % dist_trans = data_trans(:,ids1) - data_trans(:, ids2);
% % dist_trans = sqrt(ones(1,3)*(dist_trans.*dist_trans));
% % flags = double(dist_rot < Para.grid_size_rot & dist_trans < Para.grid_size_trans);
% % flags = reshape(flags, [num, num]);
% % flags = max(flags, flags');
% % [s,t,~] = find(flags);
% % ids = s < t;
% % s = s(ids);
% % t = t(ids);
% % G = graph(s,t,[],num);
% % bins = conncomp(G);
% % numClus = max(bins);
% % for clusId = 1 : numClus
% %     ids = find(bins == clusId);
% %     R_opt = zeros(3,3);
% %     t_opt = zeros(3,1);
% %     for j = 1 : length(ids)
% %         R_opt = R_opt + Result.R_opts{ids(j)};
% %         t_opt = t_opt + Result.t_opts{ids(j)};
% %     end
% %     R_opt = R_opt/length(ids);
% %     t_opt = t_opt/length(ids);
% %     [u,v,w] = svd(R_opt);
% %     R_opt = u*w';
% %     Result_clus.R_opts{clusId} = R_opt;
% %     Result_clus.t_opts{clusId} = t_opt;
% %     Result_clus.weights(clusId) = length(ids);
% % end
% % [s,order] = sort(-Result_clus.weights);
% % order = order(1:Para.K);
% % Result_clus.R_opts = Result_clus.R_opts(order);
% % Result_clus.t_opts = Result_clus.t_opts(order);
% % Result_clus.weights = Result_clus.weights(order);
% % Result_clus.weights = Result_clus.weights /sum(Result_clus.weights);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Perform clustering using a generic approach
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [Result_clus] = grid_based_clustering(Result, Para)
% % %
% % if length(Result.R_opts) <= Para.K
% %     Result_clus = Result;
% %     return;
% % end
% % num = length(Result.t_opts);
% % t_opts = zeros(3, num);
% % for id = 1 : num
% %     t_opts(:,id) = Result.t_opts{id};
% % end
% % numClus = 0;
% % clusters = translation_based_clustering(t_opts, Para);
