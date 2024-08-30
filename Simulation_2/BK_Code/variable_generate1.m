function V_cand = variable_generate1(X_train)

p = size(X_train,2);
seq = 1:1:p;
Idx = nchoosek(seq,2);
num = size(Idx,1);

V0 = ones(size(X_train,1),1);
Xl = sqrt(1.5)*(X_train-2);
Xlq = Xl;

V_cand = [V0,Xlq];

end