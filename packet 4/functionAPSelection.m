function [Krandomorder mK_AP] = functionAPSelection(M,K,Beta,nbrOfRealizations)
% Because the AP selection is done user by user in a random order, we creat
% a random vector of K user namely Krandomorder to storage the index of
% k-th user. Then we find the mk_AP suitable for k-th user by Eq (38)
% INPUT
% M         = number of the APs
% K         = number of the users
% Beta      = M x K x nbrOfRealizations large scale fading matrix
% nbrOfRealizations

%OUTPUT
% Krandomorder  = 1 x K x nbrOfRealizations matrix stores ascending vector K
% users
% mK_AP         = 1 x K x nbrOfRealizations matrix stores index of mk-th
% AP, which serves the corresponding k-th user in Krandomorder
Krandomorder    = zeros(1,K,nbrOfRealizations);
mK_AP           = zeros(1,K,nbrOfRealizations);

% Prepare for storing the random user and corresponding AP choosen in
% Krandom and mKrandom, responsively
Krandom         = zeros(1,K,nbrOfRealizations);
mKrandom         = zeros(1,K,nbrOfRealizations);

for nAPSelection = 1:nbrOfRealizations
    Krandom(:,:,nAPSelection)    = randperm(K,K);
    for kAPSelection = 1:length(Krandomorder);
        kuser = Krandom(kAPSelection);
        [maxval mkAP] = max(Beta(:,kuser,nAPSelection));
        mKrandom(1,kAPSelection,nAPSelection) = mkAP;
        Beta(mkAP,:,nAPSelection) = 0;
    end
    % Sorts the elements of Krandomorder in ascending order, then sorts mK_AP
    % responsively. Where I matrix is the sorted-value index of
    % Krandomorder matrix
    [Krandomorder I] = sort(Krandom);
    mK_AP = mKrandom(I);
end
end        