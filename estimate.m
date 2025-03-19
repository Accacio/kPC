close all;
clear all;

load('./estimation_data.mat')


[n, N]=size(x);
modes = 2^n;

% % shuffle data
% rand_perm=randperm(size(x, 2));
% x=x(:, rand_perm)
% for i=1:2
% y(:,:,i)=y(:,rand_perm,i)
% end

%%
maxIter=500;
maxErr=1e-5;
Phi0=Phibar+0.*rand(2^n,n^2+n);
% Phi0=kron([vec(A(:,:,1)).' b(:,1).'],ones(1,2^n).')+1.*rand(2^n,n^2+n);
% Phi0=[reshape(P_mult,4,4).' ]
[Phi_est,z,err,normErr]=kPC(x,y,Phi0,modes,maxIter,maxErr);
% [Phi_est,z,err,normErr]=kPC(x,y(:,:,1),[],modes,maxIter,maxErr);
Phi_est(abs(Phi_est)<1e-5)=0;
Phi_est=sortrows(Phi_est,[1 2]);

for i=1:modes
P_est(:,:,i)=reshape(Phi_est(i,1:n^2),2,2);
s_est(:,i)=reshape(Phi_est(i,n^2+1:end),1,n);
A(:,:,i)=reshape(Phibar(i,1:n^2),2,2).';
b(:,i)=reshape(Phibar(i,n^2+1:end),1,n);
end


disp(['Estimated Ps'])
disp(P_est)
disp(['Nominal Ps'])
disp(A)
% return
%% plot data
% for lambda_idx=1:2
% fig=figure;
% subplot(1,2,1)
% scatter3(x(1,:),x(2,:),y(lambda_idx,:),'k')
% sgtitle(['$\lambda_' num2str(lambda_idx) '$'],'interpreter','latex')
% end

%% final plot
for lambda_idx=1:2
fig=figure;
subplot(1,2,1)
z_colors={'r','b','k','g'};
for idx_z=1:4;
scatter3(x(1,z==idx_z),x(2,z==idx_z),y(lambda_idx,z==idx_z),z_colors{idx_z})
hold on
end
hold off
end

