using PhysicalConstants
# define the forces
function F_ab(ma,mb, r)
    r_hat = r/norm(r)
    ma_hat = ma/norm(ma)
    mb_hat = mb/norm(mb)
    return 3*μ0 * (ma * mb) /(4π * r^4)*(r_hat .* (ma_hat * mb_hat) + ma .* (r_hat * mb) + mb*(r_hat * ma))- 5 .*r_hat.*(r_hat * ma_hat) * (r_hat * mb_hat) * (ma * mb)
end
F_ba() = -F_ab()

# define the torque
