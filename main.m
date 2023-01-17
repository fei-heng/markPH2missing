%% simulated data example

clear all

tic()

% load the data set
load("sim_data.mat");

Z = [Zm, Zc];
n=size(X,1);
beta_n=size(Z,2);
nstrt=length(unique(strata));               % number of stratum

L_n=5;     % L: number of nearest neighbors
rep_n = [5]; % number of imputations

h = 0.03; % bandwidth

% parameters for tests and calculating CIF rates
a=0.01;     b=0.07;       aa=0.026;
tt=1;    zpd2=8;

% other settings
Nbindex = 2;  % Neighbor 1: H=(T,Z,A), 2: H=(T,Z)
Az_index = 2; % 1: with Az; 2: without Az

n0=17;              % Number of grid points
v=linspace(0,0.08,n0); % Grid points
K=@(x) 0.75*(1-x.^2).*(abs(x)<=1); % Kernel function


Vpred_all=[];
Vobs_all=[];
V_obs_all=[];

sel_all=[];
wei_all=[];
Z_c_all=[];


for jj=1:nstrt
    strata_sub=strata;
    strata_sub(strata~=jj-1,:)=NaN;
    X_sub=X(strata==jj-1,:);
    Z_sub=Z(strata==jj-1,:);
    Zc_sub=Zc(strata==jj-1,:);
    Zm_sub=Zm(strata==jj-1,:);
    V_sub=V(strata==jj-1,:);
    delta_sub=delta(strata==jj-1,:);
    Zaux_sub=Zaux(strata==jj-1,:);
    Vaux_sub=Vaux(strata==jj-1,:);
    
    n_jj=sum(~isnan(strata_sub));
    
    sel=~isnan(Zm_sub);
    
    
    if (Az_index == 1)
        Z_c=[Zc_sub delta_sub Zaux_sub Zc_sub.*delta_sub Zaux_sub.*delta_sub];
    else
        Z_c=[Zc_sub delta_sub Zc_sub.*delta_sub];
    end
    
    wei=sel;
    [wei_delta0]=weight(Zc_sub(delta_sub==0),sel(delta_sub==0),'parametric');
    wei(delta_sub==0)=wei_delta0;
    
    
    Vobs=V_sub;
    V_obs=V_sub(~isnan(V_sub),:);
    
    if (Nbindex == 2)
        Neigh=[ X_sub,Z_sub ];
    else
        Neigh=[ Vaux_sub, X_sub,Z_sub ];
    end
    Neigh_delta=Neigh(delta_sub==1,:);
    
    Vx=zscore(Neigh_delta);
    Vx_mis=Vx(isnan(V_sub(delta_sub==1)),:);
    Vx_obs=Vx(~isnan(V_sub(delta_sub==1)),:);
    ID=knnsearch(Vx_obs,Vx_mis,'K',L_n );
    ID_rep=[];
    for kk=1:size(Vx_mis,1)
        ID_rep(kk,:)=randsample(ID(kk,:),rep_n,true);
    end
    
    if ( size(ID_rep,1)>0 )
        Vpred0=V_obs(ID_rep);
    else
        Vpred0=[];
    end
    Vpred=repmat(Vobs,1,rep_n);
    Vpred(isnan(Vobs) & delta_sub==1,:)=Vpred0;
    
    
    
    
    Vpred_all=[Vpred_all;Vpred];
    Vobs_all=[Vobs_all;Vobs];
    V_obs_all=[V_obs_all;V_obs];
    
    sel_all=[sel_all;sel];
    wei_all=[wei_all;wei];
    Z_c_all=[Z_c_all;Z_c];
    
end


Vpred=Vpred_all;
Vobs=Vobs_all;
V_obs=V_obs_all;

sel=sel_all;
wei=wei_all;
Z_c=Z_c_all;


beta_ini=zeros(beta_n,size(v,2));%  %initial value
[eff_beta_hat,sig_eff_beta,p_values_t1,p_values_t2,...
    F10, F50, F90]=...
    aipw_imp1_strt(beta_ini,X,Z,Zc,Zm,Vpred,v,delta,sel,beta_n,...
    wei,h,n0,n,K,Z_c,rep_n,a,b,aa,strata,tt,zpd2);

F10n=reshape(F10',1,[]);
F50n=reshape(F50',1,[]);
F90n=reshape(F90',1,[]);

p_value_a=p_values_t1(2,:);
p_value_m=p_values_t1(4,:);

p_value_t2_a=p_values_t2(2,:);
p_value_t2_m=p_values_t2(4,:);

toc()

% running time: about for sample size 3000 in Alienware m15 with Intel(R) Core(TM) i7-9750H CPU @ 2.60GHz 2.59 GHz





