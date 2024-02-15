function tribool = bitcompare(a, b)
    %
    % Comparing 2 logical array in binary manner.
    %

    tribool = uint8(0);

    for idx = 1:prod(size(a))
        c = a(idx) - b(idx);

        if c
            tribool = c > 0;
            break
        end

    end

end
