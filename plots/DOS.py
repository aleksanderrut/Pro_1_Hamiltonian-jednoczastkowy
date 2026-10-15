from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt

# Ścieżka do pliku
plik = Path(
    r"C:\Users\aleks\Desktop\PWR\Numeryczne metody badania układów kwantowych\Pro_1_Hamiltonian_jednocząstkowy\data\DOS_Gauss_3D_sigma0.2_PBCfalse_Lx20_Ly20_Lz20.data"
)

# Folder do zapisu wykresów
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

match_sigma = re.search(r"_sigma([^_]+)", nazwa)
sigma = float(match_sigma.group(1)) if match_sigma else None

print("Wymiar:", wymiar)
print("PBC:", PBC)
print("Lx:", Lx)
print("Ly:", Ly)
print("Lz:", Lz)
print("Sigma:", sigma)

# WCZYTANIE DOS
data = np.loadtxt(plik)

energia = data[:, 0]
DOS = data[:, 1]

# Sortowanie po energii
indeksy = np.argsort(energia)
energia = energia[indeksy]
DOS = DOS[indeksy]

# UZUPEŁNIENIE ZERAMI DO [-8, 8]
E_min = -8.0
E_max = 8.0

# Zera po lewej stronie widma
if energia[0] > E_min:

    energia_lewa = np.array([E_min, np.nextafter(energia[0], -np.inf)])
    DOS_lewa = np.zeros(len(energia_lewa))
    energia = np.concatenate((energia_lewa, energia))
    DOS = np.concatenate((DOS_lewa, DOS))

# Zera po prawej stronie widma
if energia[-1] < E_max:

    energia_prawa = np.array([np.nextafter(energia[-1], np.inf), E_max])
    DOS_prawa = np.zeros(len(energia_prawa))
    energia = np.concatenate((energia, energia_prawa))
    DOS = np.concatenate((DOS, DOS_prawa))



# OPIS DO TYTUŁU
opis_rozmiaru = f"Lx = {Lx}"

if Ly is not None:
    opis_rozmiaru += f", Ly = {Ly}"

if Lz is not None:
    opis_rozmiaru += f", Lz = {Lz}"

if sigma is not None:
    opis_rozmiaru += rf", $\sigma$ = {sigma}"


# WYKRES
plt.figure(figsize=(9, 6))
plt.plot(energia,DOS,linestyle="-",linewidth=1.5)
plt.xlim(-8, 8)
plt.xlabel("E")
plt.ylabel("DOS(E)")
plt.title(
    f"DOS Gauss {wymiar}, PBC = {PBC}, {opis_rozmiaru}"
)
plt.grid()
plt.tight_layout()

# ZAPIS WYKRESU
plik_wykresu = folder_wykresy / f"{nazwa}.png"
plt.savefig(plik_wykresu, dpi=300)
print("Zapisano wykres:")
print(plik_wykresu)
plt.show()
