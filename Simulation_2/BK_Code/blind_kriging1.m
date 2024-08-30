function RMSE = blind_kriging1(X, M, Y, Mt, Yt)
SE = [];
for j = 1:20
    sq = 1:1:20;
    sq(j) = [];
    seq_train = sq; seq_test = j;

    X_train = kron(X(seq_train,:),ones(20,1));
    X_test = kron(X(seq_test,:),ones(51,1));
    M_train = M(:,seq_train); Y_train = Y(:,seq_train);
    M_test = Mt(:,seq_test); Y_test = Yt(:,seq_test);
    
    X_train = [X_train,M_train(:)]; Y_train = Y_train(:);
    X_test = [X_test,M_test(:)]; Y_test = Y_test(:);
    
    % scaling input variables into [1,3]
    X_train(:,1:10) = X_train(:,1:10)/250+2; X_test(:,1:10) = X_test(:,1:10)/250+2;
    X_train(:,11) = X_train(:,11)*2+1; X_test(:,11) = X_test(:,11)*2+1;
    
    % initial regression function
    V = ones(size(X_train,1),1);
    
    V_cand = variable_generate1(X_train);
    
    IDX = [12,2,11]; % chosen by Bayesian variable selection

    V = [V,V_cand(:,IDX)];
    
    model = kriging_fit1(X_train, Y_train, V);
    
    Sigy_inv = model.sigy_inv;
    
    % prediction
    Vt = ones(size(X_test,1),1);
    V_candt = variable_generate1(X_test);
    Vt = [Vt,V_candt(:,IDX)];
    
    Sig_yy = corrGaussian(model.theta(:), X_test,X_train);
    Ey = Vt*model.beta+Sig_yy*Sigy_inv*(Y_train-V*model.beta);
    se = (Ey-Y_test).^2;
    SE = [SE;se'];
end

RMSE = sqrt(mean(SE));
end