using Random
using DelimitedFiles
using LinearAlgebra

# Hamiltoanin 1D jednnocząstkowy
function gen_H_1D(L, t, PBC::Bool)
    N = L
    H = zeros(Float64, N, N)

    for i in 1:L
        if i < L
            H[i, i + 1] = -t
            H[i + 1, i] = -t
        end
    end

    if PBC == 1 && L > 2
        H[1, L] = -t
        H[L, 1] = -t
    end

    return H
end

# zwraca numerr wezła w przypadku 2D
function site_number_2D(x, y, Lx)
    return x + (y - 1) * Lx
end

# Hamiltoanin 2D jednnocząstkowy
function gen_H_2D(Lx, Ly, t, PBC::Bool)
    N = Lx * Ly
    H = zeros(Float64, N, N)
    for y in 1:Ly
        for x in 1:Lx
            i = site_number_2D(x, y, Lx)

        # przeskok w kierunku x
        if x < Lx
            j = site_number_2D(x + 1, y, Lx)
            H[i, j] = -t
            H[j, i] = -t
        elseif PBC == 1 && Lx > 2
            j = site_number_2D(1, y, Lx)
            H[i, j] = -t
            H[j, i] = -t
        end

        # przeskok w kierunku y
        if y < Ly
            j = site_number_2D(x, y + 1, Lx)
            H[i, j] = -t
            H[j, i] = -t
        elseif PBC == 1 && Ly > 2
            j = site_number_2D(x, 1, Lx)
            H[i, j] = -t
            H[j, i] = -t
        end
        end
    end
    return H
end

# zwraca numerr wezła w przypadku 3D
function site_number_3D(x, y, z, Lx, Ly)
    return x + (y - 1) * Lx + (z - 1) * Lx * Ly
end

# Hamiltoanin 3D jednnocząstkowy
function gen_H_3D(Lx, Ly, Lz, t, PBC::Bool)
    N = Lx * Ly * Lz
    H = zeros(Float64, N, N)

    for z in 1:Lz
    for y in 1:Ly
        for x in 1:Lx
            i = site_number_3D(x, y, z, Lx, Ly)

            # przeskol w kierunku x
            if x < Lx
                j = site_number_3D(x + 1, y, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            elseif PBC == 1 && Lx > 2
                j = site_number_3D(1, y, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            end

            # przeskok w kierunku y
            if y < Ly
                j = site_number_3D(x, y + 1, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t

            elseif PBC == 1 && Ly > 2
                j = site_number_3D(x, 1, z, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            end

            # przeskok w kierunku z
            if z < Lz
                j = site_number_3D(x, y, z + 1, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t
            elseif PBC == 1 && Lz > 2
                j = site_number_3D(x, y, 1, Lx, Ly)
                H[i, j] = -t
                H[j, i] = -t

            end
        end
    end
    end
    return H
end

# Zapis do pliku .data
function save_data(data, filename, PBC; seed=nothing, L_x=nothing, L_y=nothing, L_z=nothing, N_koszykow=nothing)
    folder = raw"C:\Users\aleks\Desktop\PWR\Numeryczne metody badania układów kwantowych\Pro_1_Hamiltonian_jednocząstkowy\data"
    filename_full = "$(filename)_PBC$(PBC)"
    if L_x !== nothing
        filename_full *= "_Lx$(L_x)"
    end
    if L_y !== nothing
        filename_full *= "_Ly$(L_y)"
    end
    if L_z !== nothing
        filename_full *= "_Lz$(L_z)"
    end
    if N_koszykow !== nothing
        filename_full *= "_Nkoszykow$(N_koszykow)"
    end
    if seed !== nothing
        filename_full *= "_seed$(seed)"
    end

    filename_full *= ".data"
    filepath = joinpath(folder, filename_full)

    writedlm(filepath, data)
    println("Zapisano do pliku:")
    println(filepath)
end

####################################################################################################################
function main()

    println("KOD START")

    # Tworzenie Hamiltonianu
    L_x = 90
    L_y = 90
    L_z = 20
    t = 1.0
    PBC = true

    # H1 = gen_H_1D(L_x, t, PBC)
    H2 = gen_H_2D(L_x, L_y, t, PBC)
    # H3 = gen_H_3D(L_x, L_y, L_z, t, PBC)
    
    # save_data(H1, "Hamiltonian_1D", PBC; L_x)

    ED = eigen(Symmetric(H2))
    eigenvalues = ED.values
    # eigenvectors = ED.vectors
    save_data(eigenvalues, "eigenvalues_2D", PBC; L_x, L_y)
    # save_data(eigenvectors, "eigenvectors_2D", W, PBC; seed, L_x, L_y)

    # DENSITY OF STATES - DOS
    N_koszykow = 100
    E_min = minimum(eigenvalues)
    E_max = maximum(eigenvalues)
    granice_koszykow = range(E_min, E_max, length=N_koszykow + 1)
    dE = granice_koszykow[2] - granice_koszykow[1]
    zliczenia = zeros(Int, N_koszykow)

    # robienie histogramu
    for E in eigenvalues
        przesuniecie = E - E_min
        numer_koszyka = floor(Int, przesuniecie / dE) + 1
        if numer_koszyka == N_koszykow + 1
            numer_koszyka = N_koszykow
        end
        zliczenia[numer_koszyka] += 1
    end

    # środki koszyków
    energies_DOS = zeros(N_koszykow)
    for i in 1:N_koszykow
        lewa_granica = granice_koszykow[i]
        prawa_granica = granice_koszykow[i + 1]
        energies_DOS[i] = (lewa_granica + prawa_granica) / 2
    end

    println("Liczba wartości własnych: ", length(eigenvalues))
    println("Zliczono wartości własnych: ", sum(zliczenia))

    # Zapis przed normalizacja
    histogram_data = hcat(energies_DOS, zliczenia)
    save_data(histogram_data, "Histogram_2D", PBC; L_x, L_y, N_koszykow)

    # normalizacja: suma DOS_i * dE = 1
    DOS = zliczenia ./ (length(eigenvalues) * dE)
    println("Integral DOS = ", sum(DOS) * dE) #sprawdzenie czy DOS równy 1

    # Zapis po normalizacji
    DOS_data = hcat(energies_DOS, DOS)
    save_data(DOS_data, "DOS_2D", PBC; L_x, L_y, N_koszykow)

    println("KOD ZAKOŃCZONY")

end

main()