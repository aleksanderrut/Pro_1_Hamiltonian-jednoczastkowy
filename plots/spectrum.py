from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt

# ŚCIEŻKI
plik = Path(
    r"C:\Users\aleks\Desktop\PWR\Numeryczne metody badania układów kwantowych\Pro_1_Hamiltonian_jednocząstkowy\data\eigenvalues_2D_PBCtrue_Lx90_Ly90.data"
)

folder_wykresy = Path(
    r"C:\Users\aleks\Desktop\PWR\Numeryczne metody badania układów kwantowych\Pro_1_Hamiltonian_jednocząstkowy\plots"
)

folder_wykresy.mkdir(parents=True, exist_ok=True)

# ODCZYT PARAMETRÓW Z NAZWY
nazwa = plik.stem

wymiar = re.search(r"_(\dD)_", nazwa).group(1)
PBC = re.search(r"_PBC(true|false)", nazwa).group(1)
Lx = int(re.search(r"_Lx(\d+)", nazwa).group(1))

match_Ly = re.search(r"_Ly(\d+)", nazwa)
Ly = int(match_Ly.group(1)) if match_Ly else None

match_Lz = re.search(r"_Lz(\d+)", nazwa)
Lz = int(match_Lz.group(1)) if match_Lz else None

print("Wymiar:", wymiar)
print("PBC:", PBC)
print("Lx:", Lx)
print("Ly:", Ly)
print("Lz:", Lz)

# WCZYTANIE WARTOŚCI WŁASNYCH
eigenvalues = np.loadtxt(plik)

# Sortowanie wartości własnych
eigenvalues = np.sort(eigenvalues)

# Numery wartości własnych
indeksy = np.arange(1, len(eigenvalues) + 1)

print("Liczba wartości własnych:", len(eigenvalues))
print("Minimalna energia:", eigenvalues[0])
print("Maksymalna energia:", eigenvalues[-1])

# OPIS DO TYTUŁU
opis_rozmiaru = f"Lx = {Lx}"

if Ly is not None:
    opis_rozmiaru += f", Ly = {Ly}"

if Lz is not None:
    opis_rozmiaru += f", Lz = {Lz}"

# WYKRES
plt.figure(figsize=(9, 6))
plt.plot(indeksy, eigenvalues,marker="o", markersize=2, linestyle="None")
plt.xlabel("Numer wartości własnej n")
plt.ylabel(r"Energia $E_n$")
plt.title(
    f"Wartości własne {wymiar}, PBC = {PBC}, {opis_rozmiaru}"
)
plt.grid()
plt.tight_layout()

# ZAPIS WYKRESU
plik_wykresu = folder_wykresy / f"{nazwa}.png"
plt.savefig(plik_wykresu, dpi=300)
print("Zapisano wykres:")
print(plik_wykresu)
plt.show()
