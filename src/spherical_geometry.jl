using LinearAlgebra

# NOTE: in all of these functions, the "sky plane" is the x-y plane (unit normal vector aligned with the z-axis), where relevant.

"""
Draw a random unit vector isotropically, in cartesian coordinates.
"""
function draw_random_normal_vector()
    u = 2*rand() - 1 # uniform in [-1,1]
    θ = 2π*rand()

    x = sqrt(1-u^2)*cos(θ)
    y = sqrt(1-u^2)*sin(θ)
    return [x,y,u]
end

Rx(θ) = [1 0 0; 0 cos(θ) -sin(θ); 0 sin(θ) cos(θ)] # rotation matrix around x-axis
Ry(θ) = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)] # rotation matrix around y-axis
Rz(θ) = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1] # rotation matrix around z-axis

"""
Calculate the rotation matrix that would align vector A to B.
Does not work if A = -B.
"""
function calc_rotation_matrix_A_to_B(A::Vector{Float64}, B::Vector{Float64})
    @assert(length(A) == length(B) == 3)

    v = cross(A, B)
    s = norm(v)
    c = dot(A, B)
    vx = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

    R = I + vx + ((1-c)/s^2)*vx^2
    return R
end

#calc_angle_between_vectors(A::Vector{Float64}, B::Vector{Float64}) = acos(dot(A,B)/(norm(A)*norm(B))) # angle between two vectors; WARNING: can fail if the angle is too small (near parallel vectors)!
calc_angle_between_vectors(A::Vector{Float64}, B::Vector{Float64}) = 2*atan(norm(A/norm(A) - B/norm(B)) / norm(A/norm(A) + B/norm(B))) # more stable method for computing angle between two vectors; from pg. 15 of https://people.eecs.berkeley.edu/%7Ewkahan/MathH110/Cross.pdf

"""
Calculate the orientation of a vector in a given plane (where 'nplane' is the unit normal vector to the plane).
"""
function calc_Ω_in_plane(h::Vector{Float64}, nplane::Vector{Float64})
    n = cross(nplane, h)
    Ω = acos(n[1]/norm(n))
    n[2] >= 0 ? Ω : 2π - Ω
end

# Calculate the orientation of a vector in the x-y plane (normal to z-axis).
calc_Ω_in_sky_plane(h::Vector{Float64}) = calc_Ω_in_plane(h, [0.,0.,1.])

"""
Calculate the angle between two orbits using the spherical law of Cosines.
"""
function calc_incl_spherical_cosine_law(i1::Float64, i2::Float64, ΔΩ::Float64)
    cos_i = cos(i1)*cos(i2) + sin(i1)*sin(i2)*cos(ΔΩ)
    if cos_i > 1 || cos_i < -1
        #@warn "cos_i = $(cos_i); rounding result." maxlog=10
        cos_i = round(cos_i)
    end
    return acos(cos_i)
end





"""
    calc_orbit_vector_given_system_vector(i_m, Ω, vec_ref)

Calculate the orbit normal vector given an inclination `i_m` and argument of ascending node `Ω` relative to a system reference plane, in cartesian coordinates.

# Arguments:
- `i_m::Float64`: inclination of orbit (rad) relative to the reference plane.
- `Ω::Float64`: argument of ascending node (rad) relative to the reference plane.
- `vec_ref::Vector{Float64}`: unit normal vector for the reference plane.

# Returns:
- `vec_orb::Vector{Float64}`: unit normal vector for the orbit.
"""
function calc_orbit_vector_given_system_vector(i_m::Float64, Ω::Float64, vec_ref::Vector{Float64})
    @assert(length(vec_ref) == 3)

    # Calculate rotation matrix that rotates z-axis to vec_ref:
    vec_z = [0.,0.,1.]
    inv_R = vec_z==vec_ref ? I : calc_rotation_matrix_A_to_B(vec_z, vec_ref)

    # Rotate vec_z to get the orbit normal vector:
    vec_orb = Rx(i_m)*vec_z; # rotate by i_m around x-axis
    vec_orb = Rz(Ω)*vec_orb; # rotate by Ω around z-axis
    vec_orb = inv_R*vec_orb; # rotate by inverse rotation matrix

    return vec_orb
end

"""
    calc_sky_incl_Ω_orbits_given_system_vector(i_m_list, Ω_list, vec_ref)

Calculate the inclination and orientation relative to the sky plane for all orbits in the system, given a list of mutual inclinations and arguments of ascending nodes relative to the system reference plane.

# Arguments:
- `i_m_list::Vector{Float64}`: list of orbit inclinations (rad) relative to the reference plane.
- `Ω_list::Vector{Float64}`: list of arguments of ascending node (rad) relative to the reference plane.
- `vec_ref::Vector{Float64}`: unit normal vector for the reference plane.

# Returns:
- `i_sky_list::Vector{Float64}`: list of orbit inclinations (rad) relative to the sky plane.
- `Ω_sky_list::Vector{Float64}`: list of arguments of ascending node (rad) relative to the sky plane.
"""
function calc_sky_incl_Ω_orbits_given_system_vector(i_m_list::Vector{Float64}, Ω_list::Vector{Float64}, vec_ref::Vector{Float64})
    n_pl = length(i_m_list)
    vec_z = [0.,0.,1.]

    i_sky_list = Vector{Float64}(undef, n_pl)
    Ω_sky_list = Vector{Float64}(undef, n_pl)
    for i in 1:n_pl
        vec_orb = calc_orbit_vector_given_system_vector(i_m_list[i], Ω_list[i], vec_ref)
        i_sky_list[i] = calc_angle_between_vectors(vec_z, vec_orb)
        Ω_sky_list[i] = calc_Ω_in_sky_plane(vec_orb)
    end
    return i_sky_list, Ω_sky_list
end

"""
    calc_incl_Ω_relative_to_system_invariable_plane(period_list, mass_list, ecc_list, incl_list, Ω_list; star_mass)

Calculate the inclination and orientation relative to the system invariable plane for all the orbits in the system, given a planetary system with orbital elements specified in the sky plane.

# Arguments:
- `period_list::Vector{Float64}`: a list of orbital periods (days).
- `mass_list::Vector{Float64}`: a list of planet masses (any units).
- `ecc_list::Vector{Float64}`: a list of orbital eccentricities.
- `incl_list::Vector{Float64}`: a list of inclinations (radians) relative to the sky plane.
- `Ω_list::Vector{Float64}`: a list of arguments of ascending nodes (radians) in the sky plane.
- `star_mass::Float64 = 1.`: the stellar mass (solar masses).

# Returns:
- `incl_invariable_list::Vector{Float64}`: list of orbit inclinations (rad) relative to the system invariable plane.
- `Ω_invariable_list::Vector{Float64}`: list of arguments of ascending node (rad) relative to the system invariable plane.
"""
function calc_incl_Ω_relative_to_system_invariable_plane(period_list::Vector{Float64}, mass_list::Vector{Float64}, ecc_list::Vector{Float64}, incl_list::Vector{Float64}, Ω_list::Vector{Float64}; star_mass::Float64=1.)
    n = length(period_list)
    @assert(n == length(mass_list) == length(ecc_list) == length(incl_list) == length(Ω_list) > 1)

    vec_orb_list = map(i -> calc_orbit_vector_given_system_vector(incl_list[i], Ω_list[i], [0.,0.,1.]), 1:n) # unit normals of each planet's orbital plane
    a_list = map(i -> semimajor_axis(period_list[i], star_mass), 1:n) # semi-major axes, in AU
    b_list = map(i -> a_list[i]*sqrt((1 - ecc_list[i])*(1 + ecc_list[i])), 1:n) # semi-minor axes, in AU
    L_list = map(i -> mass_list[i]*b_list[i]*sqrt(ExoplanetsSysSim.G_mass_sun_in_mks*star_mass / a_list[i]), 1:n) # angular momentum (magnitude) of each planet's orbit, as calculated from the Vis-viva equation
    Lvec_sys = sum(L_list .* vec_orb_list) # angular momentum vector of the system
    vec_invariable = Lvec_sys ./ norm(Lvec_sys) # unit normal to system invariable plane

    incl_invariable_list = map(vec_orb -> calc_angle_between_vectors(vec_invariable, vec_orb), vec_orb_list) # mutual inclinations relative to system invariable plane
    Ω_invariable_list = map(vec_orb -> calc_Ω_in_plane(vec_orb, vec_invariable), vec_orb_list) # ascending nodes relative to system invariable plane
    return incl_invariable_list, Ω_invariable_list
end


"""
    calc_incl_Ω_relative_to_system_invariable_plane(sys)

Calculate the inclination and orientation relative to the system invariable plane for all the orbits in the system, given a planetary system with orbital elements specified in the sky plane.

# Arguments:
- `sys::PlanetarySystem`: a planetary system containing the planets' physical properties (`sys.planet`) and orbital elements (`sys.orbit`) relative to the sky plane.

# Returns:
- `incl_invariable_list::Vector{Float64}`: list of orbit inclinations (rad) relative to the system invariable plane.
- `Ω_invariable_list::Vector{Float64}`: list of arguments of ascending node (rad) relative to the system invariable plane.
NOTE: the outputs default to zeros for both the inclination and ascending node if there is only one planet in the system.
"""
function calc_incl_Ω_relative_to_system_invariable_plane(sys::PlanetarySystem)
    n = length(sys.planet)
    @assert(n > 0)
    if n == 1
        incl_invariable_list, Ω_invariable_list = [0.], [0.]
    else
        period_list = map(i -> sys.orbit[i].P, 1:n)
        mass_list = map(i -> sys.planet[i].mass, 1:n)
        ecc_list = map(i -> sys.orbit[i].ecc, 1:n)
        incl_list = map(i -> sys.orbit[i].incl, 1:n)
        Ω_list = map(i -> sys.orbit[i].asc_node, 1:n)
        incl_invariable_list, Ω_invariable_list = calc_incl_Ω_relative_to_system_invariable_plane(period_list, mass_list, ecc_list, incl_list, Ω_list; star_mass=sys.star.mass)
    end
    return incl_invariable_list, Ω_invariable_list
end
