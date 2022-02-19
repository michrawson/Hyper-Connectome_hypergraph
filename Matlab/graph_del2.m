function L = graph_del2(Adj)

N1 = size(Adj,1);
N2 = size(Adj,2);

L = zeros(N1,N2);

for i = 1:N1
    L(i,:) = (sum(Adj(i,:)) + ones(1,N1)*Adj - 2*Adj(i,:))/(N1+N2-2) ...
                - Adj(i,:);
end
