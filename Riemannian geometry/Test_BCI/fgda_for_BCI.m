
%% mean class 1
clear

eps = 0.00000001; %%%

s = [16,18,20,52,55,56,89,92,93];

channel_num = length(s);

cov1 = zeros(channel_num,channel_num,0);
cov2 = zeros(channel_num,channel_num,0);

sampleRate = 1000; % Hz
lowEnd = 10; % Hz
highEnd = 30; % Hz
filterOrder = 5; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
[b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); % Generate filter coefficients

for i = 1:5
    filename = ['data_set_IVa' num2str(i) '.mat'];
    load(filename)

    data = cnt(:,s);
    setmrk = size(cnt);
    mrk.pos = [mrk.pos,setmrk(1)];
    covstru = zeros(channel_num,channel_num,280);
    for j = 1:1:(length(mrk.pos)-1)
        trail = data(mrk.pos(:,j):mrk.pos(:,j+1),:);
        trail = double(trail);
        trail = filtfilt(b, a, trail); % length*channles
        covmat = cov(trail);
        covstru(:,:,j) = covmat;%% caculate all cov matrices
    end

    class1mark = find(mrk.y==1);
    covstru1 = covstru(:,:,class1mark);
    cov1 = cat(3,cov1,covstru1);
    class2mark = find(mrk.y==2);
    covstru2 = covstru(:,:,class2mark);
    cov2 = cat(3,cov2,covstru2);
end


size1 = size(cov1);
size2 = size(cov2);

rat = 0.8; %% proportion of training samples

cov1train_num = reshape(randsample(size1(3),round(size1(3)*rat)),1,[]);
cov2train_num = reshape(randsample(size2(3),round(size2(3)*rat)),1,[]);

cov1train = cov1(:,:,cov1train_num);
cov2train = cov2(:,:,cov2train_num);

cov1va_num = setdiff(1:size1(3),cov1train_num);
cov2va_num = setdiff(1:size2(3),cov2train_num);


cov1va = cov1(:,:,cov1va_num);
cov2va = cov2(:,:,cov2va_num);
cov1ind = size(cov1va);
cov2ind = size(cov2va);

cov_va = cat(3,cov1va,cov2va);


%% FGDA
datall = cat(3,cov1train,cov2train);

data1 = cov1train;
data2 = cov2train;

[~,num_de,num_mat] = size(data1);

M = compute_cheap_mean(datall);

[V,D] = eig(M);
V_1 = V^(-1);
eigv = eig(M);

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vec1 = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = data1(:,:,i);
            Tn = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vec1(:,i) = tmp(index);
end

%% class2

[~,num_de,num_mat] = size(data2);

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vec2 = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = data1(:,:,i);
            Tn = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vec2(:,i) = tmp(index);
end



%% LDA

K = 5;
W=LDA(vec1.',vec2.',K);


%% filtering train1data
datava = cov1train;
[~,num_de,num_mat] = size(datava);


index = reshape(triu(ones(num_de)),num_de*num_de,1) == 1;

vecva = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = datava(:,:,i);
            Tn = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vecva(:,i) = tmp(index);
end

vec_new = W*inv(W.'*W)*W.'*vecva;

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vecback = zeros(1,length(index));
mat = zeros(num_de,num_de,num_mat);

w = zeros(num_de,num_de);

        for l = 1:1:num_de
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                w(l,r) = va;
            end
        end

for a = 1:num_mat
    n = 0;
    for i = 1:length(index)
        if index(i) == 1
            n = n + 1;
            vecback(i) = vec_new(n,a);
        else
             vecback(i) = 0;
        end
    end
    
    tmp = reshape(vecback,[],num_de);
    tmp2 = triu(tmp,0) + tril(tmp.',-1);
    tmp3 = V*(w.*(V_1*tmp2*V+2*D))*D*(w.*(V_1*tmp2*V+2*D))*V_1;
    mat(:,:,a) = tmp3;
end
cov1train_filtering = mat;

%% filtering train2data
datava = cov2train;
[~,num_de,num_mat] = size(datava);


index = reshape(triu(ones(num_de)),num_de*num_de,1) == 1;

vecva = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = datava(:,:,i);
            Tn = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vecva(:,i) = tmp(index);
end

vec_new = W*inv(W.'*W)*W.'*vecva;

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vecback = zeros(1,length(index));
mat = zeros(num_de,num_de,num_mat);

w = zeros(num_de,num_de);

        for l = 1:1:num_de
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                w(l,r) = va;
            end
        end

for a = 1:num_mat
    n = 0;
    for i = 1:length(index)
        if index(i) == 1
            n = n + 1;
            vecback(i) = vec_new(n,a);
        else
             vecback(i) = 0;
        end
    end
    
    tmp = reshape(vecback,[],num_de);
    tmp2 = triu(tmp,0) + tril(tmp.',-1);
    tmp3 = V*(w.*(V_1*tmp2*V+2*D))*D*(w.*(V_1*tmp2*V+2*D))*V_1;
    mat(:,:,a) = tmp3;
end

cov2train_filtering = mat;


%%  filtering validation data
datava = cov_va;
[~,num_de,num_mat] = size(cov_va);


index = reshape(triu(ones(num_de)),num_de*num_de,1) == 1;

vecva = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = datava(:,:,i);
            Tn = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vecva(:,i) = tmp(index);
end

vec_new = W*inv(W.'*W)*W.'*vecva;

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vecback = zeros(1,length(index));
mat = zeros(num_de,num_de,num_mat);

w = zeros(num_de,num_de);

        for l = 1:1:num_de
            for r = 1:1:num_de
                va =(eigv(l)+eigv(r))^(-1);
                w(l,r) = va;
            end
        end

for a = 1:num_mat
    n = 0;
    for i = 1:length(index)
        if index(i) == 1
            n = n + 1;
            vecback(i) = vec_new(n,a);
        else
             vecback(i) = 0;
        end
    end
    
    tmp = reshape(vecback,[],num_de);
    tmp2 = triu(tmp,0) + tril(tmp.',-1);
    tmp3 = V*(w.*(V_1*tmp2*V+2*D))*D*(w.*(V_1*tmp2*V+2*D))*V_1;
    mat(:,:,a) = tmp3;
end

cov_va = mat;


%%  
cov1i= ones(cov1ind(3),1);
cov2i= ones(cov2ind(3),1)*2;


covind = cat(1,cov1i,cov2i);%% index of va

class1mean = compute_bwcheap_mean(cov1train_filtering); 
class2mean = compute_bwcheap_mean(cov2train_filtering); 

%% validation

classout = zeros(length(covind),1);

for i = 1:length(covind)

   point = cov_va(:,:,i);
     dis1 = compute_W_distance(class1mean,point);
     dis2 = compute_W_distance(class2mean,point);
   if dis1 > dis2
       classout(i) = 2;
   else 
       classout(i) = 1;
   end 
end
lengt = size(covind);
accuracyrate = (lengt(1)-sum(abs(classout-covind)))/lengt(1)


