from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt


# Ścieżka do pliku
plik = Path(
    r"C:\Users\aleks\Desktop\PWR\Numeryczne metody badania układów kwantowych\Pro_1_Hamiltonian_jednocząstkowy\data\DOS_1D_PBCfalse_Lx8000_Nkoszykow200.data"
)

# Folder do zapisu wykresów
folder_wykresy = Path(
    r"C:\Users\aleks\Desktop\PWR\Numeryczne metody badania układów kwantowych\Pro_1_Hamiltonian_jednocząstkowy\plots"
)


# =========================
# ODCZYT PARAMETRÓW Z NAZWY
# =========================

nazwa = plik.stem

wymiar = re.search(r"_(\dD)_", nazwa).group(1)
PBC = re.search(r"_PBC(true|false)", nazwa).group(1)
Lx = int(re.search(r"_Lx(\d+)", nazwa).group(1))

match_Ly = re.search(r"_Ly(\d+)", nazwa)
Ly = int(match_Ly.group(1)) if match_Ly else None

match_Lz = re.search(r"_Lz(\d+)", nazwa)
Lz = int(match_Lz.group(1)) if match_Lz else None

match_Nkoszykow = re.search(r"_Nkoszykow(\d+)", nazwa)
Nkoszykow = int(match_Nkoszykow.group(1)) if match_Nkoszykow else None

print("Wymiar:", wymiar)
print("PBC:", PBC)
print("Lx:", Lx)
print("Ly:", Ly)
print("Lz:", Lz)
print("Nkoszykow:", Nkoszykow)


# =========================
# WCZYTANIE DOS
# =========================

data = np.loadtxt(plik)

energia = data[:, 0]
DOS = data[:, 1]

dE = energia[1] - energia[0]


# =========================
# UZUPEŁNIENIE ZERAMI DO [-8, 8]
# =========================

energia_lewa = np.arange(-8, energia[0], dE)
DOS_lewa = np.zeros(len(energia_lewa))

energia_prawa = np.arange(energia[-1] + dE, 8 + dE / 2, dE)
DOS_prawa = np.zeros(len(energia_prawa))

energia = np.concatenate((energia_lewa, energia, energia_prawa))
DOS = np.concatenate((DOS_lewa, DOS, DOS_prawa))


# =========================
# OPIS DO TYTUŁU
# =========================

opis_rozmiaru = f"Lx = {Lx}"
if Ly is not None:
    opis_rozmiaru += f", Ly = {Ly}"
if Lz is not None:
    opis_rozmiaru += f", Lz = {Lz}"
if Nkoszykow is not None:
    opis_rozmiaru += f", N koszyków = {Nkoszykow}"


# =========================
# WYKRES
# =========================

plt.plot(
    energia,
    DOS,
    marker="o",
    markersize=3
)

plt.xlim(-8, 8)
plt.xlabel("E")
plt.ylabel("DOS(E)")
plt.title(f"DOS {wymiar}, PBC = {PBC}, {opis_rozmiaru}")
plt.grid()
plt.tight_layout()


# =========================
# ZAPIS WYKRESU
# =========================

plik_wykresu = folder_wykresy / f"{nazwa}.png"
plt.savefig(plik_wykresu,dpi=300)
print("Zapisano wykres:")
print(plik_wykresu)
plt.show()