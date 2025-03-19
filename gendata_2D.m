% Example of estimation data generation

%= Clean variables
close all
clear

paren = @(x, varargin) x(varargin{:}); %
                                       % apply index in created elements
curly = @(x, varargin) x{varargin{:}}; %

%%
n=2;
values=linspace(0,20,10);
[ v{1:n} ]=ndgrid(values);
x(:,:) =cell2mat(cellfun(@(x) reshape(x,[],1),v,'UniformOutput',0))';
A=zeros(n,n,2^n);
b=zeros(n,1,2^n);
% TODO(accacio) create automaticly
A(:,:,1)=[1 0; 0 1];
b(:,:,1)=[2;3];
A(:,:,2)=[1 4; 6 3];
b(:,:,2)=[5;9];
A(:,:,3)=[2 0; 5 8];
b(:,:,3)=[7;3];
A(:,:,4)=[1 8; 6 9];
b(:,:,4)=[8;3];
y=zeros(n,length(x));
cond1=[1 0]*x<10;
cond2=[0 2]*x<30;
part1=~cond1&~cond2;
part2=~cond1&cond2;
part3=cond1&~cond2;
part4=cond1&cond2;
y(:,part1)=A(:,:,1)*x(:,part1)+b(:,:,1);
y(:,part2)=A(:,:,2)*x(:,part2)+b(:,:,2);
y(:,part3)=A(:,:,3)*x(:,part3)+b(:,:,3);
y(:,part4)=A(:,:,4)*x(:,part4)+b(:,:,4);

%% plot
for y_idx=1:2
fig=figure;
subplot(1,2,1)
scatter3(x(1,:),x(2,:),y(y_idx,:,1),'k')
sgtitle(['$y_' num2str(y_idx) '$'],'interpreter','latex')
end

Phibar=[vec(A(:,:,1)).' b(:,:,1).';
    vec(A(:,:,2)).' b(:,:,2).'
    vec(A(:,:,3)).' b(:,:,3).'
    vec(A(:,:,4)).' b(:,:,4).'];
Phibar=sortrows(Phibar);
save('estimation_data.mat','x','y','Phibar');
