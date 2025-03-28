function A = combineArrs(A,B)
    bIdx = find(B(1,:));
    aFirstZero = find(~A(1,:),1);
    A(:, aFirstZero:aFirstZero+length(bIdx)-1) = B(:,bIdx);
end