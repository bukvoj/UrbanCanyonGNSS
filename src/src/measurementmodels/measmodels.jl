function meas(x,u,p,t)
    # p is: (vector of positions, vector of sv vels, validity flags, biases)
    svpos = p[1]
    svvel = p[2]
    valid = p[3]
    biases = p[4]
    
    numsvs = length(svpos)
    l2available = (length(biases) == 2*numsvs)
    if l2available
        y = zeros(4*numsvs)
    else
        y = zeros(2*numsvs)
    end
    
    for i in 1:numsvs
        y[i] = norm(x[1:3] - svpos[i]) + x[7]               + biases[i]
        if l2available
            y[i + numsvs] = norm(x[1:3] - svpos[i]) + x[7]      + biases[i+numsvs]
        end
    end
    for i in 1:numsvs
        los = (svpos[i]-x[1:3]) / norm(svpos[i]-x[1:3])
        y[i + 2*numsvs] = -los'*(svvel[i] - x[4:6]) - x[8]
        if l2available
            y[i + 3*numsvs] = -los'*(svvel[i] - x[4:6]) - x[8]
        end
    end
    y[isnan.(y)] .= 0
    return y .* valid
end

function meas_jac(x,u,p,t)
    # p is: (vector of positions, vector of sv vels, validity flags, biases)
    
    numsvs = length(p[1])
    svpos = p[1]
    valid = p[3]

    l2available = (length(p[4]) == 2*numsvs)
    if l2available
        J = zeros(4*numsvs,8)
        for i in 1:numsvs
            J[i,1:3] = (x[1:3] - svpos[i])' / norm(x[1:3] - svpos[i])
            J[i,7] = 1
            J[i + numsvs,1:3] = (x[1:3] - svpos[i])' / norm(x[1:3] - svpos[i])
            J[i + numsvs,7] = 1
        end
        for i in 1:numsvs
            J[2*numsvs + i,4:6] = -(x[1:3] - svpos[i])' / norm(x[1:3] - svpos[i])
            J[2*numsvs + i,8] = -1
            J[3*numsvs + i,4:6] = -(x[1:3] - svpos[i])' / norm(x[1:3] - svpos[i])
            J[3*numsvs + i,8] = -1
        end
    else
        J = zeros(2*numsvs,8)
        for i in 1:numsvs
            J[i,1:3] = (x[1:3] - svpos[i])' / norm(x[1:3] - svpos[i])
            J[i,7] = 1
        end
        for i in 1:numsvs
            J[numsvs + i,4:6] = -(x[1:3] - svpos[i])' / norm(x[1:3] - svpos[i])
            J[numsvs + i,8] = -1
        end
    end
    J[isnan.(J)] .= 0
    J .*= valid
    return J
end