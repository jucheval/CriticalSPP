function vec_to_mat_test(v::Vector)
    d = Int((-1 + isqrt(1 + 8 * length(v))) ÷ 2)
    output = zeros(d, d)

    idx = 1
    for i in 1:d
        output[i, i] = v[idx]
        idx += 1
    end

    for i in 1:(d - 1)
        for j in (i + 1):d
            output[i, j] = v[idx]
            output[j, i] = v[idx] # Symmetric entry
            idx += 1
        end
    end
    return output
end
