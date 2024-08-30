function V_cand = variable_generate(X_train)

p = size(X_train,2);
seq = 1:1:p;
Idx = nchoosek(seq,2);
num = size(Idx,1);

V0 = ones(size(X_train,1),1);
Xl = sqrt(1.5)*(X_train-2); % Xq = sqrt(4.5)*(X_train-2).^2-sqrt(2);
Xlq = Xl; % Xlq = [Xl,Xq];

V1 = [];
for i = 1:num
    idx1 = Idx(i,1);
    idx2 = Idx(i,2);
    if abs(idx2-idx1)~=p
        V = Xlq(:,idx1).*Xlq(:,idx2);
        V1 = [V1,V];
    end
end
V_cand = [V0,Xlq,V1];

end