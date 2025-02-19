function [Phi,z,err,norm_err] = kPC(X,Y,phi_init,modes,maxIter,maxErr)
% Bradley, P., & Mangasarian, O. (2000). K-Plane Clustering. J Global Optim, 16(1), 23–32. http://dx.doi.org/10.1023/a:1008324625522

    k = modes;
    [n, N] = size(Y);
    z = zeros(1,N);

    if(isempty(phi_init))
        Phi=20*rand(modes,n^2+n);
    else
        Phi=phi_init;
    end
    OldPhi=zeros(size(Phi));
    exit_msg=0;

    for j=1:maxIter
        %= kPC Cluster assignement
        Upsilon=kron(ones(n,1),eye(n));
        Delta=kron(eye(N),ones(1,n));
        G=kron(ones(1,N),eye(n));
        Y_=kron(G,ones(n,1));
        Omega=[(Upsilon*X*Delta).*Y_; G].';
        y_cell = cellfun(@(x) x',mat2cell(Y',ones(1,N))','UniformOutput',0);
        y_diag=blkdiag(y_cell{:});

        % calculate errors
        for i=1:modes
            y_lin=reshape(Omega*-Phi(i,:)',n,N);
            y_lin_cell = cellfun(@(x) x',mat2cell(y_lin',ones(1,N))','UniformOutput',0);
            y_lin_diag=blkdiag(y_lin_cell{:});
            err_diag=y_diag-y_lin_diag;
            err_cell=y_lin_cell-y_cell;

            norm_err_cell(i,:)=cellfun(@(x) norm(x),err_cell,'UniformOutput',0);
        end
        err=cell2mat(err_cell);
        norm_err=cell2mat(norm_err_cell);
        [~, z] = min(norm_err);
        Responsibilities = kron(1:modes,ones(1,N).').'==z;

        %= Cluster Update
        for i=1:modes
            responsibilities=Responsibilities(i,:);
            resp2=cellfun(@(x) x*eye(n),mat2cell(responsibilities',ones(1,N)),'UniformOutput',0);
            Gamma=sqrt(sparse(blkdiag(resp2{:})));

            Phi(i,:)=-(Gamma*Omega)\(Gamma*Y(:));
        end
        OldPhi=Phi;

        if norm(err,'fro')<maxErr
            exit_msg=1;
            break;
        end

        % if norm(OldPhi-Phi,'fro')<maxErr
        %     exit_msg=2;
        %     break;
        % end



    end

    switch exit_msg
    case 1
      disp(['maxErr residues'])
    case 2
        disp(['minimal iteration on estimated'])
    otherwise
      disp(['max iter reached'])
    end

end
