%% import edf.file
[tt,anno] = edfread('/Users/jc/Documents/MATLAB/M 20/EEG data/files/S100/S100R03.edf');

%% transform edf

hz = length(tt.Fc5_{1,1});
[timelength,num_chann]= size(tt);
T = timetable2table(tt);
tablecell= table2cell(T);
tablecell = tablecell(:,2:65);

table_str = zeros(timelength*hz,num_chann);
for j = 1:num_chann
for i = 1:timelength
    table_str((1+hz*(i-1)):(hz*i),j) = tablecell{i,j};
end
end
%%
% anno

Anno = timetable2table(anno);
atablecell= table2cell(Anno);
[num_row,num_col] = size(atablecell);

tab = table('Size',[num_row 2],'VariableTypes',{'string','double'}); 
tab(:,1) = atablecell(:,2);
for i  = 1:num_row
[x,origUnit] = time2num(atablecell{i,1},seconds);
tab.Var2(i) = x;
[y,origUnit] = time2num(atablecell{i,3},seconds);
tab.Var3(i) = y;
end

%% T1
Index1 = find(tab.Var1 == 'T1');
mat_arr_1 = zeros(64,64,length(Index1));


for i = 1:length(Index1)
    
num_r = round(tab.Var2(Index1(length(Index1)))*hz+tab.Var3(Index1(length(Index1)))*hz);

if num_r > timelength*hz
    num_r = timelength*hz;
end
    mat = table_str(round(tab.Var2(Index1(i))*hz):num_r,:);
    mat_arr_1(:,:,i) = cov(mat);

end

%% T2
Index2 = find(tab.Var1 == 'T2');
mat_arr_2 = zeros(64,64,length(Index2));


for i = 1:length(Index2)
    
num_r = round(tab.Var2(Index2(length(Index2)))*hz+tab.Var3(Index2(length(Index2)))*hz);

if num_r > timelength*hz
    num_r = timelength*hz;
end
    mat = table_str(round(tab.Var2(Index2(i))*hz):num_r,:);
    mat_arr_2(:,:,i) = cov(mat);

end
%% collect file name
fileFolder=fullfile('/Users/jc/Documents/MATLAB/M 20/EEG data/files/S001');
 
dirOutput=dir(fullfile(fileFolder,'*.edf'));
 
fileNames={dirOutput.name};

% collect filepaths

maindir='/Users/jc/Documents/MATLAB/M 20/EEG data/files/';
subdir=dir(maindir);
subdir(1:3,:) = [];
%len=length(subdir);
vect = [1,7,8,29,34,35,42,55,60,62,70,71,72,85,93];
len =length(vect);
acctable = zeros(len,8);


for num_c = vect

%num_c = 2;
% Task 1 06

 filepath=fullfile(maindir,subdir(num_c).name,strcat(subdir(num_c).name,'R06.edf'));    
 [cov_1,cov_2] = import_cov_data(filepath);
 cov_1_a = cov_1;
 cov_2_a = cov_2;

% Task 1 10


 filepath=fullfile(maindir,subdir(num_c).name,strcat(subdir(num_c).name,'R10.edf'));    
 [cov_1,cov_2] = import_cov_data(filepath);
 cov_1_b = cov_1;
 cov_2_b = cov_2;

% Task 1 14

 filepath=fullfile(maindir,subdir(num_c).name,strcat(subdir(num_c).name,'R14.edf'));    
 [cov_1,cov_2] = import_cov_data(filepath);
 cov_1_c = cov_1;
 cov_2_c = cov_2;

% Concatenate cov arrays
cov1 = cat(3,cov_1_a,cov_1_b,cov_1_c);
cov2 = cat(3,cov_2_a,cov_2_b,cov_2_c);




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


datall = cat(3,cov1train,cov2train);
data1 = cov1train;
data2 = cov2train;



% FGDA class 1
[~,num_de,num_mat] = size(data1);

M = compute_bwpro_mean(datall);

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

% class2

[~,num_de,num_mat] = size(data2);

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vec2 = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = data2(:,:,i);
            Tn = (M * p)^(1/2)+(p * M)^(1/2) - 2 * M;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vec2(:,i) = tmp(index);
end

% LDA

K = 60;
W=LDA(vec1.',vec2.',K);


% filtering train1data
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

% filtering train2data
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


%  filtering validation data
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


%  
cov1i= ones(cov1ind(3),1);
cov2i= ones(cov2ind(3),1)*2;


covind = cat(1,cov1i,cov2i);%% index of va

class1mean = compute_bwpro_mean(cov1train_filtering); 
class2mean = compute_bwpro_mean(cov2train_filtering); 

% validation

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

lengt_bwfd = size(covind);


accuracyrate_bwfd = (lengt_bwfd(1)-sum(abs(classout-covind)))/lengt_bwfd(1);

acctable(num_c,1)= lengt_bwfd(1);
acctable(num_c,2) = accuracyrate_bwfd;





%% %%%%%%%%%%%% Rie %%%%%%%%%%%%%%%
datall = cat(3,cov1train,cov2train);

data1 = cov1train;
data2 = cov2train;

[~,num_de,num_mat] = size(data1);

M = rie_mean(datall,15);


A = M ^ (1/2);    
B = A ^ (-1);    

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vec1 = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = data1(:,:,i);
            Tn = A * logm(B * p * B) * A;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vec1(:,i) = tmp(index);
end

% class2

[~,num_de,num_mat] = size(data2);

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vec2 = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = data2(:,:,i);
            Tn = A * logm(B * p * B) * A;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vec2(:,i) = tmp(index);
end



% LDA

K = 60;

W=LDA(vec1.',vec2.',K);


% filtering train1data
datava = cov1train;
[~,num_de,num_mat] = size(datava);


index = reshape(triu(ones(num_de)),num_de*num_de,1) == 1;

vecva = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = datava(:,:,i);
           Tn = A * logm(B * p * B) * A;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vecva(:,i) = tmp(index);
end

vec_new = W/(W.'*W)*W.'*vecva;

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vecback = zeros(1,length(index));
mat = zeros(num_de,num_de,num_mat);

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
    tmp3 = A * expm(B * tmp2 * B) * A; 
    mat(:,:,a) = tmp3;
end
cov1train_filtering = mat;

% filtering train2data
datava = cov2train;
[~,num_de,num_mat] = size(datava);


index = reshape(triu(ones(num_de)),num_de*num_de,1) == 1;

vecva = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = datava(:,:,i);
            Tn = A * logm(B * p * B) * A;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vecva(:,i) = tmp(index);
end

vec_new = W*inv(W.'*W)*W.'*vecva;

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vecback = zeros(1,length(index));
mat = zeros(num_de,num_de,num_mat);


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
    tmp3 = A * expm(B * tmp2 * B) * A; 
    mat(:,:,a) = tmp3;
end

cov2train_filtering = mat;


%  filtering validation data
datava = cov_va;
[~,num_de,num_mat] = size(cov_va);


index = reshape(triu(ones(num_de)),num_de*num_de,1) == 1;

vecva = zeros(num_de*(num_de-1)/2+num_de,num_mat);

for i = 1:1:num_mat
            p = datava(:,:,i);
            Tn = A * logm(B * p * B) * A;
            tmp = reshape(triu(Tn),num_de*num_de,1);
            vecva(:,i) = tmp(index);
end

vec_new = W*inv(W.'*W)*W.'*vecva;

index = reshape(triu(ones(num_de)),num_de*num_de,1)==1;

vecback = zeros(1,length(index));
mat = zeros(num_de,num_de,num_mat);


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
    tmp3 = A * expm(B * tmp2 * B) * A; 
    mat(:,:,a) = tmp3;
end

cov_va = mat;



covind = cat(1,cov1i,cov2i);%% index of va

class1mean = compute_riepro_mean(cov1train_filtering,'A'); 
class2mean = compute_riepro_mean(cov2train_filtering,'A'); 

% validation

classout = zeros(length(covind),1);

for i = 1:length(covind)

   point = cov_va(:,:,i);
     dis1 = rie_dis(class1mean,point);
     dis2 = rie_dis(class2mean,point);
   if dis1 > dis2
       classout(i) = 2;
   else 
       classout(i) = 1;
   end 
end

lengt_rfd = size(covind);
accuracyrate_rfd = (lengt_rfd(1)-sum(abs(classout-covind)))/lengt_rfd(1);

acctable(num_c,3)= lengt_rfd(1);
acctable(num_c,4) = accuracyrate_rfd;

%% Directly  validation Bw pro

cov1i= ones(cov1ind(3),1);
cov2i= ones(cov2ind(3),1)*2;

covind = cat(1,cov1i,cov2i);%% index of va

class1mean = compute_bwpro_mean(data1); 
class2mean = compute_bwpro_mean(data2); 

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
lengt_b = size(covind);
accuracyrate_b = (lengt_b(1)-sum(abs(classout-covind)))/lengt_b(1);
acctable(num_c,5)= lengt_b(1);
acctable(num_c,6) = accuracyrate_b;
%% %% Directly  validation Rie pro

cov1i= ones(cov1ind(3),1);
cov2i= ones(cov2ind(3),1)*2;

covind = cat(1,cov1i,cov2i);%% index of va

class1mean = compute_riepro_mean(data1,'A'); 
class2mean = compute_riepro_mean(data2,'A'); 

classout = zeros(length(covind),1);

for i = 1:length(covind)

   point = cov_va(:,:,i);
     dis1 = rie_dis(class1mean,point);
     dis2 = rie_dis(class2mean,point);
   if dis1 > dis2
       classout(i) = 2;
   else 
       classout(i) = 1;
   end 
end
lengt_r = size(covind);
accuracyrate_r = (lengt_r(1)-sum(abs(classout-covind)))/lengt_r(1);
acctable(num_c,7)= lengt_r(1);
acctable(num_c,8) = accuracyrate_r;

end

%%
acc_bf = sum(acctable(:,1).*acctable(:,2))/sum(acctable(:,1))
acc_rf = sum(acctable(:,3).*acctable(:,4))/sum(acctable(:,3))
acc_b = sum(acctable(:,5).*acctable(:,6))/sum(acctable(:,5))
acc_r = sum(acctable(:,7).*acctable(:,8))/sum(acctable(:,7))

%% git example 1

% Data formating
COVtest = data.data(:,:,data.idxTest);
trueYtest  = data.labels(data.idxTest);

COVtrain = data.data(:,:,data.idxTraining);
Ytrain  = data.labels(data.idxTraining);

%% MDM classification - Multiclass
metric_mean = {'euclid','logeuclid','riemann','ld'};
metric_dist = {'euclid','logeuclid','riemann','ld','kullback'};
acc = zeros(length(metric_mean),length(metric_dist));

for i=1:length(metric_mean)
    for j=1:length(metric_dist)
        Ytest = mdm(COVtest,COVtrain,Ytrain,metric_mean{i},metric_dist{j});
        acc(i,j) = 100*mean(Ytest==trueYtest);
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows : distance metric, Colums : mean metric');
disp('------------------------------------------------------------------');
displaytable(acc',metric_mean,10,{'.1f'},metric_dist)
disp('------------------------------------------------------------------');

%% Discriminant geodesic filtering + MDM classification - Multiclass
metric_mean = {'euclid','logeuclid','riemann','ld'};
metric_dist = {'euclid','logeuclid','riemann','ld','kullback'};
acc = zeros(length(metric_mean),length(metric_dist));

for i=1:length(metric_mean)
    for j=1:length(metric_dist)
        Ytest = fgmdm(COVtest,COVtrain,Ytrain,metric_mean{i},metric_dist{j});
        acc(i,j) = 100*mean(Ytest==trueYtest);
    end
end
disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows : distance metric, Colums : mean metric');
disp('------------------------------------------------------------------');
displaytable(acc',metric_mean,10,{'.1f'},metric_dist)
disp('------------------------------------------------------------------');

%% MDM classification - Binary case
metric_mean = 'riemann';
metric_dist = 'riemann';
acc = diag(nan(4,1));

for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = mdm(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),Ytrain(ixtrain),metric_mean,metric_dist);
        % Accuracy
        acc(i,j) = 100*mean(Ytest==trueYtest(ixtest));
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%% Discriminant geodesic filtering + MDM Classification - Binary case
metric_mean = 'riemann';
metric_dist = 'riemann';
acc = diag(nan(4,1));

for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = fgmdm(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),Ytrain(ixtrain),metric_mean,metric_dist);
        % Accuracy
        acc(i,j) = 100*mean(Ytest==trueYtest(ixtest));
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%% Kmeans usupervised Classification - Binary case
metric_mean = 'riemann';
metric_dist = 'riemann';
acc = diag(nan(4,1));

% for each couple of classes
for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = kmeanscov(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),2,metric_mean,metric_dist);
        % Find the right labels
        Classes = unique(trueYtest(ixtest));
        truelabels = (trueYtest(ixtest) == Classes(1))+1;
        % Accuracy
        acc(i,j) = 100*mean(Ytest==truelabels);
        if acc(i,j)<50
            acc(i,j) = 100-acc(i,j);
        end
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%% Tangent Space LDA Classification - Binary case
% the riemannian metric
metric_mean = 'riemann';
% update tangent space for the test data - necessary if test data corresponds to
% another session. by default 0.
update = 0;
acc = diag(nan(4,1));

for i=1:4
    for j=i+1:4
        % Select the trials
        ixtrain = (Ytrain==i)|(Ytrain==j);
        ixtest = (trueYtest==i)|(trueYtest==j);
        % Classification
        Ytest = tslda(COVtest(:,:,ixtest),COVtrain(:,:,ixtrain),Ytrain(ixtrain),metric_mean,update);
        % Accuracy
        acc(i,j) = 100*mean(Ytest==trueYtest(ixtest));
    end
end

disp('------------------------------------------------------------------');
disp('Accuracy (%) - Rows/Colums : Couple of classes');
disp('------------------------------------------------------------------');
displaytable(acc'+acc,{'Right Hand','Left Hand','Foot','Tongue'},10,{'.1f'},{'Right Hand','Left Hand','Foot','Tongue'})
disp('------------------------------------------------------------------');

%%
 [le,ns]= size(mat);

sampleRate = 128; % Hz
lowEnd = 10; % Hz
highEnd = 30; % Hz
filterOrder = 5; % Filter order (e.g., 2 for a second-order Butterworth filter). Try other values too
[b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2));
mat = double(mat);
trail = filter(b, a, mat);


Fs = 128;

y0 = abs(fft(trail(:,2))); 
f = (0:le-1)*Fs/le;
plot(f,y0);
xlabel('Frequency'); 
ylabel('Amplitude');
%%
y = bandpass(mat,[10 30],Fs);
y0 = abs(fft(y(:,1))); 
plot(f,y0);
xlabel('Frequency'); 
ylabel('Amplitude');
%%
vect = [1,7,8,29,34,35,42,55,60,62,70,71,72,85,93];
