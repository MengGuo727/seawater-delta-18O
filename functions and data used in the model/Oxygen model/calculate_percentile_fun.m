function[data_5,data_25,data_50,data_75,data_95] = ...
    calculate_percentile_fun(data,nt,t,size)
% this function returns the middle 50%, 90%, and the median 
% values of the given dataset
%
% Meng Guo, Yale University
% Summer, 2019

% create empty matrixes to store the results
data_5 = zeros(1,size);
data_25 = zeros(1,size);
data_50 = zeros(1,size);
data_75 = zeros(1,size);
data_95 = zeros(1,size);

% sort the given dataset
for i= 1:nt
    data_sort(i,:) = sort(data(i,:));
end

% get the middle 50%, 90%, and median values of the given dataset
for i = 1:nt
    data_25(i) = data_sort(i,max(1,round(size*0.25)));
    data_75(i) = data_sort(i,max(1,round(size*0.75)));
    data_5(i) = data_sort(i,max(1,round(size*0.05)));
    data_95(i) = data_sort(i,max(1,round(size*0.95)));
    data_50(i) = data_sort(i,max(1,round(size*0.50)));
end