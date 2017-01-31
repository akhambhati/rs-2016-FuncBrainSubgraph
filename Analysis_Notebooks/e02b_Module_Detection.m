%% Read Data and Format for Processing
clear; clc;

load ~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.AdjMatr.mat

% Convert adj_matr to cell array
for s=1:length(adj_matr)
    A{s} = squeeze(adj_matr(s, :, :));
end
clear s adj_matr;

% Convert adj_name to cell array
for s=1:length(adj_name)
    Label{s} = strsplit(adj_name(s, :), '.');
end
clear s adj_name;

%% Generate modularity matrix (categorical for func handle)
gamma = 1.0;
omega = 1.0;

N=length(A{1});
T=length(A);
ii=[]; jj=[]; vv=[];
twomu=0;
for s=1:T
    indx=[1:N]'+(s-1)*N;
    [i,j,v]=find(A{s});
    ii=[ii;indx(i)]; jj=[jj;indx(j)]; vv=[vv;v];
    k=sum(A{s});
    kv=zeros(N*T,1);
    twom=sum(k);
    twomu=twomu+twom;
    kv(indx)=k/twom;
    kcell{s}=kv;
end
AA = sparse(ii,jj,vv,N*T,N*T);
clear ii jj vv
kvec = full(sum(AA));
all2all = N*[(-T+1):-1,1:(T-1)];
AA = AA + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
twomu=twomu+T*omega*N*(T-1);
save('~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.ModMat_FnHandle.mat', ...
     'AA', 'kcell', 'kvec', 'N', 'T', 'gamma', 'omega', 'twomu', '-v7.3')

%% Perform genlouvain (categorical for func handle)
clear; clc; 

load ~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.ModMat_FnHandle.mat

B = @(i) AA(:,i) - gamma*kcell{ceil(i/(N+eps))}*kvec(i);
[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);
imagesc(S);

%% Generate modularity matrix (categorical for full matrix)
gamma = 1.0;
omega = 5.0;

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
save('~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.ModMat_Full.mat', ...
     'B', 'N', 'T', 'gamma', 'omega', 'twomu', '-v7.3')
 
 %% Perform genlouvain (categorical for full matrix)
clear; clc; 

load ~/Remotes/RSRCH.NMF_Subnetworks/e02b-DynFuncModule-Population/Module_Optimization.ModMat_Full.mat

[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);
imagesc(S);