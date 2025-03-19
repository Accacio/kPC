% Example of estimation data generation from dMPC article

%= Clean variables
close all
clear

paren = @(x, varargin) x(varargin{:}); %
                                       % apply index in created elements
curly = @(x, varargin) x{varargin{:}}; %

%% Options

%= Optimization settings
options = optimset('Display', 'off');
warning off
% options = [];
rand('seed',2);

%= Simulation parameters
M=4;    %= # of systems
Te=.25; %= Sampling
Np=2;   %= Prediction horizon

simK = 1;    %= Simulation horizon
negotP = 200; %= max # of iteration for each negotiation
err_theta=1e-3; %= err to test theta convergence in each negotiation

% TODO(accacio): change name of these variables
chSetpoint_list = 1; %= Change setpoint?
selfish_list = 1;    %= Selfish behaviour?
secure_list = 1;     %= Secure algorithm?
%= Global constraint
% Umax=7;
Umax=7;

%= Input bounds
u_min=0;
% u_min=-inf;
u_max=Umax;
% u_max=inf;

%% Define systems
% TODO(accacio): move to a function

Cair_mean=8;
Cwalls_mean=5;
Roaia_mean=5;
Riwia_mean=2.5;
Rowoa_mean=1.1;

Cair  = repmat(Cair_mean,1,M)+(-.5+rand(1,M));
Cwalls = repmat(Cwalls_mean,1,M)+(-.5+rand(1,M));
Roaia = repmat(Roaia_mean,1,M)+(-.5+rand(1,M));
Riwia = repmat(Riwia_mean,1,M)+(-.5+rand(1,M));
Rowoa = repmat(Rowoa_mean,1,M)+(-.5+rand(1,M));

%= Define continuos systems using 3R2C
for i=M:-1:1 % make it backward to "preallocate"
    csys(:,:,1,i)=model_3R2C(Roaia(i),Riwia(i),Rowoa(i),Cwalls(i),Cair(i));
end
ni=size(csys.B(:,:,1,1),2); %= # of inputs
ns=size(csys.A(:,:,1,1),2); %= # of states
n=Np*ni; % constant used everywhere
dsys=c2d(csys,Te); %= Discretize systems

%= search space
% see https://accacio.gitlab.io/blog/matlab_combinations/
values=linspace(0,2*Umax,10);
[ v{1:n} ]=ndgrid(values);
em_theta(:,:) =cell2mat(cellfun(@(x) reshape(x,[],1),v,'UniformOutput',0))';


%= Output prediction matrices, such $\vec{Y}=\mathcal{M}\vec{x}_k+\mathcal{C}\vec{U}$
% These functions use operator overloading, see https://accacio.gitlab.io/blog/matlab_overload/
Mmat_fun=@(sys,n) ...
         cell2mat( ...
    (repmat(mat2cell(sys.C,1),1,n) .* ...
     (repmat(mat2cell(sys.A,size(sys.A,2)),1,n).^num2cell(1:n)) ...
    )'  ...
                 );
Cmat_fun=@(sys, n) ...
         cell2mat( ...
    paren( ...
    (repmat(mat2cell(sys.C, 1), 1, n+1) .* ...
     (horzcat( ...
    zeros(size(sys.A)), ...
    repmat(mat2cell(sys.A, size(sys.A,2)),1,n).^num2cell(1:n)...
             )) .* ...
     repmat(mat2cell(sys.B, size(sys.B,1), size(sys.B,2)), 1, n+1)) ...
    , tril(toeplitz(1:n))+1));

%= H and f, such $\frac{1}{2}\vec{U}^TH\vec{U}+f^T\vec{U}$
H_fun=@(Cmat,Q,R) round(Cmat'*Q*Cmat+R*eye(size(Q)),10);
f_fun=@(Cmat,Mmat,Q,xt,Wt) Cmat'*Q*(Mmat*xt-Wt);
c_fun=@(Mmat,Q,Xt,Wt) Xt'*Mmat'*Q*Mmat*Xt-2*Wt'*Q*Mmat*Xt+Wt'*Q*Wt;

%= Gains Q and R for $\sum_{j=1}^n \|v\|^2_{Q}+\|u\|^2_{R}$
for i=M:-1:1 % make it backward to "preallocate"
    Qbar(:,:,i)=10*eye(Np*size(dsys(:,:,1,i).C,1)); % no x no
    Rbar(:,:,i)=eye(n); % nc x nc
end

%= Prediction matrices for the systems
for i=M:-1:1
    Mmat(:,:,i)=Mmat_fun(dsys(:,:,1,i),Np);
    Cmat(:,:,i)=Cmat_fun(dsys(:,:,1,i),Np);
    H(:,:,i)=H_fun(Cmat(:,:,i),Qbar(:,:,i),Rbar(:,:,i));
end

Ac = kron(ones(M,1),eye(n)); %
                             % Global Constraints
bc = kron(ones(Np,1),Umax);  %

% clear -regexp [^f].*_fun % Delete all functions but f_fun
clear Cwalls* Cair* R* csys

%= Initial state
% TODO(accacio): move up in code
% TODO(accacio): use random variables?
X0(:,1,1) = [17 3.2]';
X0(:,1,2) = [20 5.3]';
X0(:,1,3) = [15 3.1]';
X0(:,1,4) = [17. 5.7]';

%= Setpoint
Wt = [X0(1,1)*1.20;
      X0(1,2)*1.20;
      X0(1,3)*1.20;
      X0(1,4)*1.20;
     ];
Wt_final = [Wt(1)*1.05;
            Wt(2)*1.05;
            Wt(3)*1.05;
            Wt(4)*1.05;
           ];
Wt_change_time =[ simK/2;
                  simK;
                  simK;
                  simK;
                ];

umin(1:M)=u_min;
umax(1:M)=u_max;

%= Selfish behavior profile (cheating matrix)
T(:,:,1)=1*eye(n);
% T(:,:,1)=1*diag(rand(n,1));
T(:,:,2)=1*eye(n);
T(:,:,3)=1*eye(n);
T(:,:,4)=1*eye(n);

%= Time selfish behavior activated
selfish_time= [ simK/2;
                simK;
                simK;
                simK;
              ];

rho_fun = @(a,b,negot) 1/(a+b*negot);
% rho = @(a,b,negot) (1/a)/negot;
a=100;
b=100;

%% === Control Loop ===
for chSetpoint=chSetpoint_list
for selfish=selfish_list
for secure=secure_list
tic

u=zeros(n,M);
J=zeros(simK,M);
theta=zeros(n,negotP,simK,M);
lambda=zeros(n,M);
lambdaHist=zeros(n,negotP,simK,M);
uHist=zeros(ni,simK,M);
xt=zeros(ns,simK,M);
xt(:,1,1:M)=X0(:,1,1:M);
lastp=zeros(simK);
norm_err=zeros(simK,M);
PhiHist=zeros(2^n,n^2+n,simK,M);

for k=1:simK
    k
    rho=rho_fun(a,b,k);
    %= update setpoint?
    % TODO(accacio): add chSetpoint
    if (chSetpoint)
        for i=M:-1:1
            if(k>Wt_change_time(i))
                Wt(i)=Wt_final(i);
            end
        end
    end

    %= Get value of f[k]
    for i=M:-1:1
        f(:,i)=f_fun(Cmat(:,:,i),Mmat(:,:,i),Qbar(:,:,i),xt(:,k,i),Wt(i));
        fHist(:,k,i)=f(:,i);
        cHist(k,i)=c_fun(Mmat(:,:,i),Qbar(:,:,i),xt(:,k,i),kron(ones(Np,1),Wt(i)));
    end

    % %= EM


    em_lambda=zeros(n,size(em_theta,2),M);
    em_u=zeros(n,M);
    em_J=zeros(n,M);
    if secure
    for i=1:M
        for cur_theta=1:size(em_theta,2)
            % QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0)
            [~,~,~,~,l] = quadprog(H(:,:,i), f(:,i), ...
                                              eye(ni*n), em_theta(:,cur_theta), ...
                                              [], [], ...
                                              umin(:,i)*ones(ni*n,1), ...  % Lower Bound
                                              umax(:,i)*ones(ni*n,1), ...  % Upper Bound
                                              [], options);
            em_lambda(:,cur_theta,i)=l.ineqlin;
            if selfish
                if k>selfish_time(i)
                    em_lambda_mod(:,cur_theta,i)=T(:,:,i)*em_lambda(:,cur_theta,i);
                end
            end
        end
    end

    end

end
toc
end
end
end
%%
x=em_theta;
y=em_lambda(:,:,1);

P = -H(:,:,1);
s = -f(:,1);
invP=inv(P);

for k=1:2^n-1 % 0-th is already empty
active=nonzeros(bitget(k,1:n).*(1:n));
inactive=setdiff(1:n,active);
P_mult(active,active,k+1)=adjoint(invP(active,active)).'./det(P(inactive,inactive))*det(P);
% BUG(accacio) correct s_mult
s_mult(:,k+1)=f(:,1);
end

Phibar=[vec(P_mult(:,:,1)).' f(:,1).';
    vec(P_mult(:,:,2)).' f(:,2).'
    vec(P_mult(:,:,3)).' f(:,3).'
    vec(P_mult(:,:,4)).' f(:,4).'];
Phibar=sortrows(Phibar);
save('estimation_data.mat','x','y','Phibar');
