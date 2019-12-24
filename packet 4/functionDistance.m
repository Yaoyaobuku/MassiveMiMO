function [d_MK xM yM xK yK] = functionDistance(M, K, D_sqr, DistanceControl, nbrOfRealizations)
% Creat the (M x 1 x nbrOfRealizations) matrix for storing the xy-coordinate of M-APs 
xM = zeros(M, 1, nbrOfRealizations);
yM = zeros(M, 1, nbrOfRealizations);
    
% Creat the (K x 1 x nbrOfRealizations) matrix for storing the
% xy-coordinate of K-users
xK = zeros(K, 1, nbrOfRealizations);
yK = zeros(K, 1, nbrOfRealizations);

% Creat the (M x K x nbrOfRealizations) matrix for storing the distance
% between m-th AP and k-th user
d_MK = zeros(M, K, nbrOfRealizations);

for nDistance = 1:nbrOfRealizations
    if (strcmp(DistanceControl, 'Halton') == 1)
        %% Solution A
        % Create a Halton sequence point set
        P = haltonset(2); 
        %constructs a new Halton sequence point set object in 2-dimensions. First dimension 
        % presents x-coordinate; and y-coordinate is showed by the remain piece.
    
        % Creat (M + K) uniformly distributed points
        XY = D_sqr*net(P, M+K); 
        % (M + K by 2) matrix denotes xy-coordinates of (M + K) points
        % where D_sqrt is the size of the rectangular covering (this code only use D_sqr because D =
        % 1000m x 1000m. Generally, for the A x B rectangular we should
        % multiple A or B by the corresponding vector of P.
    
        % Extract the coordinates of M APs and K users
        %     if ((length(XY(:,1)) == length(unique(XY(:,1))))&(length(XY(:,2)) == length(unique(XY(:,2))))== 1)
        % Check the same coordinate points. Actually, the variables in Halton
        % sequence are unique, therefore don't need checking the same coordinate
        % points step.
        xM(:, 1, nDistance) = XY(1:M, 1);
        yM(:, 1, nDistance) = XY(1:M, 2);
        xK(:, 1, nDistance) = XY(M+1:M+K, 1);
        yK(:, 1, nDistance) = XY(M+1:M+K, 2);
    else
        %% Solution B: Using Uniform Distribution
        % Creat (M + K) x-coordinates of (M + K) AP-s and users
        pdX = makedist('Uniform',0,D_sqr); % where [0 D_sqr] is the coordinate range
        pdY = makedist('Uniform',0,D_sqr);
    
        % Extract the coordinates of M APs and K users
        X = random(pdX, M+K, 1);
        Y = random(pdY, M+K, 1);
    
        xM(:, 1, nDistance) = X(1:M, 1);
        yM(:, 1, nDistance) = Y(1:M, 1);
        xK(:, 1, nDistance) = X(M+1:M+K, 1);
        yK(:, 1, nDistance) = Y(M+1:M+K, 1);
    end
    
    % Calculate the distance between m-th AP and k-th user at nDistance
    % iteration
    for mDistance = 1:M
        for kDistance = 1:K
            d_MK(mDistance, kDistance, nDistance) = sqrt([xM(mDistance, 1, nDistance) ...
             - xK(kDistance, 1, nDistance)]^2 + [yM(mDistance, 1, nDistance) - yK(kDistance, 1, nDistance)]^2);
        end
    end
end