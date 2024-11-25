function [Xr,NN_MF]=MC_NeuMF1(X,M,s,options,k)
%%% X:[n m]
% Code for DBDNMF: A Dual Branch Deep Neural Matrix Factorization Method for Drug Response Prediction
% inputs
% -- X: an n by m matrix with missing entries (m variables, n samples)
% -- M: binary mask matrix for X, 1 indicate observed, 0 indicate missed
% -- s: network structure, a row vector with L elements, indicating one
%       input layer, L-2 hidden layers, and one output layer;
%       e.g. [r 5*r 10*r m], the last value must equal to 'm'; r<<m
% -- options:
%    options.Wp: the weight decay penalty, e.g. 1e-3
%    options.Zp: the latent variable penalty, e.g. 1e-3
%    options.activation_func: the activation functions for each layer
%       the element should be one of 'tanh_opt','softmax', 'sigm', 'linear'.
%       the length should be L-1, corresponding to s.
%       e.g., {'tanh_opt','linear'}.
%    k: the weight of nonlinear structure in the network.
% outputs
% -- Xr: the recovered matrix
% -- NN_MF: the network structure and parameters

[n,m]=size(X);
d=s(1);
NN_MF=NN_MF_setup(s,options);
NN_MF.s=s;
NN_MF.Z=zeros(n,d);
NN_MF.m=m;
NN_MF.n=n;
NN_MF.d=d;
NN_MF.Wp=options.Wp;
NN_MF.Zp=options.Zp;
NN_MF.MC=1;
NN_MF.M=M;
NN_MF.X=X;
NN_MF.k=k;
disp('Training neural networks for DMFMC ......')
NN_MF.maxiter=options.maxiter;
NN_MF=NN_MF_optimization_rprop(NN_MF);
Xr=X.*M+(k*NN_MF.a{4}+(1-k)*NN_MF.a{5}).*~M;
NN_MF.Z=NN_MF.a{1}(:,2:end);
%
end

%% optimize by RPROP 
function NN_MF=NN_MF_optimization_rprop(NN_MF)
p.verbosity = 3;                    % Increase verbosity to print something [0~3]
p.MaxIter   = NN_MF.maxiter;          
p.d_Obj     = 1e-5;
p.method    = 'IRprop+';          
p.display   = 0;

w=[];
for i=1:length(NN_MF.W)
    w=[w;NN_MF.W{i}(:)];
end
z=NN_MF.Z(:);
y=[z;w];
[y,f,EXITFLAG,STATS] = rprop(@fg_DMFMC,y,p,NN_MF);
lz=NN_MF.d*NN_MF.n;
Z=reshape(y(1:lz),NN_MF.n,NN_MF.d);
t=lz+1;
for i=1:length(NN_MF.W)
    [a,b]=size(NN_MF.W{i});
    NN_MF.W{i}=reshape(y(t:t+a*b-1),a,b);
    t=t+a*b;
end
NN_MF=NN_MF_ff(NN_MF,Z);
%
end


%% f, dW, dZ
function [f,g]=fg_DMFMC(y,NN_MF)
    lz=NN_MF.d*NN_MF.n;
    Z=reshape(y(1:lz),NN_MF.n,NN_MF.d);
    t=lz+1;
    for i=1:length(NN_MF.W)
        [a,b]=size(NN_MF.W{i});
        NN_MF.W{i}=reshape(y(t:t+a*b-1),a,b);
        t=t+a*b;
    end
    NN_MF = NN_MF_ff(NN_MF,Z);
    NN_MF = NN_MF_bp(NN_MF);
    gW=[];
    w=y(lz+1:end);
    sum_w=0;
    for i=1:length(NN_MF.W)
        wt=NN_MF.W{i}(:,2:end);
        sum_w=sum_w+sum(wt(:).^2);
        dW=NN_MF.dW{i}+NN_MF.Wp*[zeros(size(NN_MF.W{i},1),1) NN_MF.W{i}(:,2:end)];
        gW=[gW;dW(:)];
%         gW=[gW;NN_MF.dW{i}(:)+myNN.weight_penalty_L2*myNN.W{i}(:)];
    end
    gZ=NN_MF.dZ+NN_MF.Zp*Z/NN_MF.n;
    gZ=gZ(:);
    f=NN_MF.loss+0.5*NN_MF.Wp*sum_w+0.5*NN_MF.Zp*sum(Z(:).^2)/NN_MF.n;
    g=[gZ;gW];

    %W_all = (NN_MF.W{3}.*NN_MF.W{4})*[zeros(1,size(NN_MF.W{1},2));(NN_MF.W{1}.*NN_MF.W{2})];
    %writematrix(W_all,'W_all.csv')
    

end


%% 
function NN_MF=NN_MF_setup(s,options)
if length(options.activation_func)~=2*(length(s)-1)
    error('The number of layers does not match the number of activation functions!')
end
NN_MF.activation_func=options.activation_func;
NN_MF.layer=length(s)-1;

NN_MF.W{1}=(rand(s(2),s(1)+1)-0.5)*2*4*sqrt(6/(s(2)+s(1)));
NN_MF.W{2}=(rand(s(2),s(1)+1)-0.5)*2*4*sqrt(6/(s(2)+s(1)));

NN_MF.W{3}=(rand(s(3),s(2)+1)-0.5)*2*4*sqrt(6/(s(3)+s(2)));
NN_MF.W{4}=(rand(s(3),s(2)+1)-0.5)*2*4*sqrt(6/(s(3)+s(2)));
%for i=1:NN_MF.layer
%    NN_MF.W{i}=(rand(s(i+1),s(i)+1)-0.5)*2*4*sqrt(6/(s(i+1)+s(i)));
%end
end
%%
function NN_MF=NN_MF_ff(NN_MF,Z)
Z=[ones(NN_MF.n,1) Z];% add bias 1
L=length(NN_MF.s);
NN_MF.a{1}=Z; 

switch NN_MF.activation_func{1}
    case 'sigm'
        NN_MF.a{2}=sigm(NN_MF.a{1}*NN_MF.W{1}');
    case 'tanh_opt'
        NN_MF.a{2}=tanh_opt(NN_MF.a{1}*NN_MF.W{1}');
    case 'relu_opt'
        NN_MF.a{2}=relu_opt(NN_MF.a{1}*NN_MF.W{1}');
    case 'linear'
        NN_MF.a{2}=NN_MF.a{1}*NN_MF.W{1}';
end
NN_MF.a{2}=[ones(NN_MF.n,1) NN_MF.a{2}];

switch NN_MF.activation_func{2}
    case 'sigm'
        NN_MF.a{3}=sigm(NN_MF.a{1}*NN_MF.W{2}');
    case 'tanh_opt'
        NN_MF.a{3}=tanh_opt(NN_MF.a{1}*NN_MF.W{2}');
    case 'relu_opt'
        NN_MF.a{3}=relu_opt(NN_MF.a{1}*NN_MF.W{2}');
    case 'linear'
        NN_MF.a{3}=NN_MF.a{1}*NN_MF.W{2}';
end
NN_MF.a{3}=[ones(NN_MF.n,1) NN_MF.a{3}];

switch NN_MF.activation_func{3}
    case 'sigm'
        NN_MF.a{4}=sigm(NN_MF.a{2}*NN_MF.W{3}');
    case 'tanh_opt'
        NN_MF.a{4}=tanh_opt(NN_MF.a{2}*NN_MF.W{3}');
    case 'relu_opt'
        NN_MF.a{4}=relu_opt(NN_MF.a{2}*NN_MF.W{3}');
    case 'linear'
        NN_MF.a{4}=NN_MF.a{2}*NN_MF.W{3}';
end

switch NN_MF.activation_func{4}
    case 'sigm'
        NN_MF.a{5}=sigm(NN_MF.a{3}*NN_MF.W{4}');
    case 'tanh_opt'
        NN_MF.a{5}=tanh_opt(NN_MF.a{3}*NN_MF.W{4}');
    case 'relu_opt'
        NN_MF.a{5}=relu_opt(NN_MF.a{3}*NN_MF.W{4}');
    case 'linear'
        NN_MF.a{5}=NN_MF.a{3}*NN_MF.W{4}';
end

% pedictive error and value of loss function
NN_MF.e=NN_MF.X-(NN_MF.k*NN_MF.a{4}+(1-NN_MF.k)*NN_MF.a{5});
if NN_MF.MC==1
    NN_MF.e=NN_MF.e.*NN_MF.M;
end
switch NN_MF.activation_func{4}
    case {'sigm', 'linear','tanh_opt','relu_opt'}
        NN_MF.loss=1/2*sum(sum(NN_MF.e.^2))/NN_MF.n;
%     case 'softmax'
%         NN_MF.loss=-sum(sum(Y.* log(NN_MF.a{L})))/NN_MF.n;
end
%
end
%%
function NN_MF=NN_MF_bp(NN_MF)
L = length(NN_MF.s);

switch NN_MF.activation_func{3}
    case 'sigm'
        d{4}=-NN_MF.e.*(NN_MF.a{4}.*(1-NN_MF.a{4}));
    case'relu_opt'
        d{4}=-NN_MF.e.*(NN_MF.a{4} > 0);
    case 'tanh_opt'
        d{4}=-NN_MF.e.*(1.7159*2/3 *(1-1/(1.7159)^2*NN_MF.a{4}.^2));
    case {'softmax','linear'}
        d{4}=-NN_MF.e;
end

switch NN_MF.activation_func{4}
    case 'sigm'
        d{5}=-NN_MF.e.*(NN_MF.a{5}.*(1-NN_MF.a{5}));
    case'relu_opt'
        d{5}=-NN_MF.e.*(NN_MF.a{5} > 0);
    case 'tanh_opt'
        d{5}=-NN_MF.e.*(1.7159*2/3 *(1-1/(1.7159)^2*NN_MF.a{5}.^2));
    case {'softmax','linear'}
        d{5}=-NN_MF.e;
end

switch NN_MF.activation_func{1}
   case 'sigm'
       d_act1=NN_MF.a{2}.*(1-NN_MF.a{2});
   case'relu_opt'
       d_act1=NN_MF.a{2} > 0;
   case 'linear'
       d_act1=1;
   case 'tanh_opt'
       d_act1=1.7159*2/3 *(1-1/(1.7159)^2*NN_MF.a{2}.^2);
end

switch NN_MF.activation_func{2}
   case 'sigm'
       d_act2=NN_MF.a{3}.*(1-NN_MF.a{3});
   case'relu_opt'
       d_act1=NN_MF.a{3} > 0;
   case 'linear'
       d_act2=1;
   case 'tanh_opt'
       d_act2=1.7159*2/3 *(1-1/(1.7159)^2*NN_MF.a{3}.^2);
end

d{2} = (d{4} * NN_MF.W{3}).* d_act1;
d{3} = (d{5} * NN_MF.W{4}).* d_act2;% Bishop (5.56)

%

NN_MF.dW{3} = (d{4}' * NN_MF.a{2}) / size(d{4}, 1);
NN_MF.dW{4} = (d{5}' * NN_MF.a{3}) / size(d{5}, 1);

NN_MF.dW{1} = (d{2}(:,2:end)' * NN_MF.a{1}) / size(d{2}, 1);      
NN_MF.dW{2} = (d{3}(:,2:end)' * NN_MF.a{1}) / size(d{3}, 1);

% dZ
if size(d{2},2)>size(NN_MF.W{1},1)
    NN_MF.dZ=NN_MF.k*d{2}(:,2:end)*NN_MF.W{1}(:,2:end)+...
             (1-NN_MF.k)*d{3}(:,2:end)*NN_MF.W{2}(:,2:end);
else
    NN_MF.dZ=NN_MF.k*d{2}(:,1:end)*NN_MF.W{1}(:,2:end)+...
             (1-NN_MF.k)*d{3}(:,1:end)*NN_MF.W{2}(:,2:end);
end
NN_MF.dZ=NN_MF.dZ/NN_MF.n;
end