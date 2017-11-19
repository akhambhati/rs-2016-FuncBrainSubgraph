%% Read Data and Format for Processing
clear; clc;

load ~/Remotes/RSRCH.NMF_Subnetworks/e02b1-DynFuncModule-Population/Module_Optimization.AdjMatr.mat
%adj_matr = adj_matr(1:80*5,:,:);

% Convert adj_matr to cell array
for s=1:length(adj_matr)
    A{s} = squeeze(adj_matr(s, :, :));
end

% Convert adj_name to cell array
for s=1:length(adj_matr)
    Label{s} = strsplit(adj_name(s, :), '.');
end

% Find new_subject indices
new_subject(1) = 0;
for s=2:length(adj_matr)
    if Label{s-1}{1} == Label{s}{1}
        new_subject(s) = 0;
    else
        new_subject(s) = 1;
    end
end
time_subj = find(new_subject)-1;
clear s adj_matr adj_name;

%% Generate modularity matrix (categorical for full matrix)
gamma = 1.0;
omega = 1.0;
rho = 0.00001;

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);

% Link time points
Tsubj = [1, N*time_subj, N*T];
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
for ii=2:length(Tsubj)-1
    B(Tsubj(ii)+1:Tsubj(ii)+1+N, Tsubj(ii)+1-N:Tsubj(ii)+1) = 0;
    B(Tsubj(ii)+1-N:Tsubj(ii)+1, Tsubj(ii)+1:Tsubj(ii)+1+N) = 0;
end

% Link subjects
all2all = N*[(-T+1):-1,1:(T-1)];
B_all2all = rho*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
for ii=1:length(Tsubj)-1
    B_all2all(Tsubj(ii):Tsubj(ii+1), Tsubj(ii):Tsubj(ii+1)) = 0;    
end

% Add back subjects
B = B + B_all2all;

save('~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.ModMat_Full.mat', ...
     'B', 'N', 'T', 'time_subj', 'gamma', 'omega', 'twomu', '-v7.3')
 
 %% Perform genlouvain (categorical for full matrix)
clear; clc;

load ~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.ModMat_Full.mat

[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S, N, T);

time_subj = [0, time_subj, T];
for ii=1:length(time_subj)-1
    Ssubj{ii} = S(:, time_subj(ii)+1:time_subj(ii+1));
end    

for ii=1:length(Ssubj)
    figure();
    imagesc(Ssubj{ii});
end