clear
load('MPRC_SZ_JimWaltz_6_8_20.mat');

size(H)
size(P_minusone)

for i = 1:124
    cor_NC(:, :, i) = corr(H(:, :, i));
end;

for i = 1:103
    cor_SZ(:, :, i) = corr(P_minusone(:, :, i));
end;

size(cor_NC)
size(cor_SZ)

cor_SZ_fisherZ = double.empty();
cor_SZ_fisherZ_permute = double.empty();
for i = 1:103
    cor_SZ_fisherZ(:, :, i) = atanh(cor_SZ(:, :, i));
    cor_SZ_fisherZ_temp = cor_SZ_fisherZ(:, :, i);
    cor_SZ_fisherZ_temp(logical(eye(size(cor_SZ_fisherZ_temp)))) = 0;
    cor_SZ_fisherZ_permute(i, :) = squareform(cor_SZ_fisherZ_temp);
end;

cor_NC_fisherZ = double.empty();
cor_NC_fisherZ_permute = double.empty();
for i = 1:124
    cor_NC_fisherZ(:, :, i) = atanh(cor_NC(:, :, i));
    cor_NC_fisherZ_temp = cor_NC_fisherZ(:, :, i);
    cor_NC_fisherZ_temp(logical(eye(size(cor_NC_fisherZ_temp)))) = 0;
    cor_NC_fisherZ_permute(i, :) = squareform(cor_NC_fisherZ_temp);
end;

Q30_SZ = quantile(cor_SZ_fisherZ_permute, 0.3, 2);
Q30_NC = quantile(cor_NC_fisherZ_permute, 0.3, 2);
Q30_all = vertcat(Q30_SZ, Q30_NC);
Q30_diff = Q30_all - mean(Q30_all(:));

cor_SZ_fisherZ_permute_correct = double.empty();
cor_NC_fisherZ_permute_correct = double.empty();
for i=1:103
    cor_SZ_fisherZ_permute_correct(i, :) = cor_SZ_fisherZ_permute(i,:) - Q30_diff(i);
end
for i=1:124
    j = i + 103;
    cor_NC_fisherZ_permute_correct(i, :) = cor_NC_fisherZ_permute(i,:) - Q30_diff(j);
end


cor_SZ_fisherZ_correct = double.empty();
cor_NC_fisherZ_correct = double.empty();
for i = 1:103
    cor_SZ_fisherZ_correct(:, :, i) = squareform(cor_SZ_fisherZ_permute_correct(i, :));
end;
for i = 1:124
    cor_NC_fisherZ_correct(:, :, i) = squareform(cor_NC_fisherZ_permute_correct(i, :));
end;

save('MPRC_transformed_cor.mat','cor_SZ_fisherZ_correct','cor_NC_fisherZ_correct')

