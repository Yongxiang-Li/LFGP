function R = corr_generate1(X_train,theta)

p = size(X_train,2);

Rl = [];
for i = 1:p
    rl = (3-3*phi(2,theta(i)))/(3+4*phi(1,theta(i))+2*phi(2,theta(i)));
    Rl = [Rl,rl];
end

Rlq = Rl;

R = diag([1,Rlq]);


    function corr = phi(x,para)
        corr = exp(-para*x^2);
    end

end