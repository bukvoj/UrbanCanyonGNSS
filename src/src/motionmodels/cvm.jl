function cvm(x,u,p,t)
    A = cvm_jac(x,u,p,t)
    return A*x
end

function cvm_jac(x,u,p,t)
    Δt = 1
    A = [1 0 0 Δt 0 0 0 0;0 1 0 0 Δt 0 0 0;0 0 1 0 0 Δt 0 0;0 0 0 1 0 0 0 0;0 0 0 0 1 0 0 0;0 0 0 0 0 1 0 0;0 0 0 0 0 0 1 Δt;0 0 0 0 0 0 0 1]
    return A
end