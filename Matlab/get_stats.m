function [accuracy, F1, Rsq_class] = get_stats(Y_pred_class, Y_class)

accuracy = 100*length(Y_pred_class(Y_pred_class==Y_class))/length(Y_pred_class);

Y_pred_class = double(Y_pred_class);
Y_class = double(Y_class);

true_pos = sum((Y_pred_class == 2) & (Y_class == 2));
false_pos = sum((Y_pred_class == 2) & (Y_class == 1));
false_neg = sum((Y_pred_class == 1) & (Y_class == 2));
true_neg = sum((Y_pred_class == 1) & (Y_class == 1));

% true_pos = sum((Y_pred_class == 1) & (Y_class == 1));
% false_pos = sum((Y_pred_class == 1) & (Y_class == 2));
% false_neg = sum((Y_pred_class == 2) & (Y_class == 1));
% true_neg = sum((Y_pred_class == 2) & (Y_class == 2));

if true_pos+false_pos > 0
    precision = true_pos/(true_pos+false_pos);
else
    precision = 0;
end

assert(isfinite(precision));

if true_pos+false_neg > 0
    recall = true_pos/(true_pos+false_neg);
else
    recall = 0;
end

assert(isfinite(recall));

if precision+recall > 0
    assert(isfinite(2*precision*recall/(precision+recall)));
    F1 = 2*precision*recall/(precision+recall);
else
    F1 = 0;
end

res_ss = sum((Y_class-Y_pred_class).^2);

tot_ss = sum((Y_class-mean(Y_class)).^2);

Rsq_class = 1-res_ss/tot_ss;
