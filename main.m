clear all;close all;clc;

% import dataset

Y_filepath ='dataset\sim1\truecounts.csv';
[Y, text, alldata] = xlsread(Y_filepath);
Y_dim = size(Y);

[mask, text, alldata] = xlsread('dataset\sim1\dropout.csv');

Known_index = find(mask==0);
Missing_index = find(mask==1);

[Y_ture, text, alldata] = xlsread('dataset\sim1\truecounts.csv');
Y_ture(Known_index)=0;

Y_0 = Y;
Y_0(Missing_index) = 0;


Mask_mat = ones(size(Y));
Mask_mat = logical(Mask_mat);

Y_recover1 = zeros(Y_dim);
result = zeros(30,9);
for cv_run = 1:1
    disp(sprintf('Run %d:',cv_run));
    rand('seed',cv_run + 20000);
    

    Y_tmp = Y_0;
    Mask_test = ~mask;
    M = (Mask_test - Mask_mat)+1;
    
    m=2000;
    n=500;
    r=500;
    
    % scDBImpute setup
    s=[r 600 m];% input size, hidden size 1, ..., output size

    options.Wp=0.01;
    options.Zp=0.01;
    options.maxiter=3500;
    % 'tanh_opt','sigm','linear'
    options.activation_func={'tanh_opt','linear','linear','linear'};

    k = 0.2;
    [Y_DMF,NN_MF]=scDBImpute(Y_tmp',Mask_test',s,options,k);

    Yr = Y_DMF';
    Yr(Yr<0)=0;

    % compute recovery error
    re_error=norm((Y_ture-Yr).*(1-M),'fro')/norm(Y_ture.*(1-M),'fro');
    
    disp(['Relative recovery error is ' num2str(re_error)]) ;
    Y_recover1(Missing_index) = Yr(Missing_index);
end


    
    figure;
    subplot(2,2,1);imagesc(Y_0);colorbar;
    subplot(2,2,2);imagesc(Y_recover1);colorbar;
    subplot(2,2,3);scatter(Y_ture(Missing_index),Y_recover1(Missing_index),'.');
    global_pcc = corr(Y_ture(Missing_index),Y_recover1(Missing_index));
    global_err = sqrt(sum((Y_ture(Missing_index)-Y_recover1(Missing_index)).^2)/length(~isnan(Y_ture(Missing_index))));
    cosSim = dot(Y_ture(Missing_index),Y_recover1(Missing_index))/(norm(Y_ture(Missing_index))*norm(Y_recover1(Missing_index)));
    L1 = pdist2(median(Y_ture(Missing_index)),median(Y_recover1(Missing_index)),'cityblock');
    

    
