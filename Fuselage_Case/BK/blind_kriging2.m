function RMSE = blind_kriging2(X, Xt, M, Y, Mt, Yt)
SE = [];
sq = 1:1:30;
seq_train = sq;

X_train = kron(X(seq_train,:),ones(60,1));
M_train = M(:,seq_train); Y_train = Y(:,seq_train);
X_train = [X_train,M_train(:)]; Y_train = Y_train(:);


% scaling input variables into [1,3]
X_train(:,1:10) = X_train(:,1:10)/450+2;
X_train(:,11) = X_train(:,11)*2+1;

% initial regression function
V = ones(size(X_train,1),1);  
V_cand = variable_generate1(X_train);

% Bayesian variable selection
% IDX = [];
% for i = 1:3
%     model = kriging_fit2(X_train, Y_train, V);
%     R = corr_generate1(X_train,model.theta);
%     Sigy_inv = model.sigy_inv;
%     u = R*V_cand'*Sigy_inv*(Y_train-V*model.beta);
%     [~,I] = max(abs(u));
%     V = [V,V_cand(:,I)];
%     IDX = [IDX,I];
% end

tic;
IDX = [12,2,11]; % chosen by Bayesian variable selection
V = [V,V_cand(:,IDX)];
model = kriging_fit2(X_train, Y_train, V);
Sigy_inv = model.sigy_inv;
t = toc;
fprintf('the running time for each prediction is %7.4f s\n',t);

for j = 1:25
    seq_test = j;
    X_test = kron(Xt(seq_test,:),ones(60,1));
    M_test = Mt(:,seq_test); Y_test = Yt(:,seq_test);
    X_test = [X_test,M_test(:)]; Y_test = Y_test(:);

    % scaling input variables into [1,3]
    X_test(:,1:10) = X_test(:,1:10)/450+2;
    X_test(:,11) = X_test(:,11)*2+1;

    % prediction
    Vt = ones(size(X_test,1),1);
    V_candt = variable_generate1(X_test);
    Vt = [Vt,V_candt(:,IDX)];
    
    Sig_yy = corrGaussian(model.theta(:), X_test,X_train);
    Ey = Vt*model.beta+Sig_yy*Sigy_inv*(Y_train-V*model.beta);
    se = (Ey-Y_test).^2;
    SE = [SE;se'];

end
RMSE = sqrt(mean(SE,2));
end