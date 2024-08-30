function R = corr_generate(X_train,theta)

p = size(X_train,2);
seq = 1:1:p;
Idx = nchoosek(seq,2);
num = size(Idx,1);

Rl = [];
for i = 1:p
    rl = (3-3*phi(2,theta(i)))/(3+4*phi(1,theta(i))+2*phi(2,theta(i)));
    Rl = [Rl,rl];
end

% Rq = [];
% for j = 1:p
%     rq = (3-4*phi(1,theta(j))+phi(2,theta(j)))/(3+4*phi(1,theta(j))+2*phi(2,theta(j)));
%     Rq = [Rq,rq];
% end

Rlq = Rl; % Rlq = [Rl,Rq];

R1 = [];
for i = 1:num
    idx1 = Idx(i,1);
    idx2 = Idx(i,2);
    if abs(idx2-idx1)~=p
        r = Rlq(idx1)*Rlq(idx2);
        R1 = [R1,r];
    end
end
R = diag([1,Rlq,R1]);


    function corr = phi(x,para)
        corr = exp(-para*x^2);
    end

end