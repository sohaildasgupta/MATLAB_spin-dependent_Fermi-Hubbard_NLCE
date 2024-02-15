function cycle = findpermutation(a, b)
    %
    % Find ordering indexing that tranforms a to b, i.e. a(cycle)=b.
    %
    %
    len = length(a);
    [~, s1] = sort(a);
    [~, s2] = sort(b);
    Id = eye(len);
    cycle = (Id(s2, :))' * Id(s1, :) * (1:len)';
end
