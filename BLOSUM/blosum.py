import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp

#/run/media/changery/0123-4567/ENCB/Bioinformatica/Proyecto/Secuencias/fasta_envelope/DENV3_MUSCLE.fas

# Funciones

def get_col(ncol: int, seqs: list):
    col = ""
    for seq in seqs:
        col += seq[ncol]
    return col


# Aminoacidos globales
amin = (
    '-',
    'A',
    'C',
    'D',
    'E',
    'F',
    'G',
    'H',
    'I',
    'K',
    'L',
    'M',
    'N',
    'P',
    'Q',
    'R',
    'S',
    'T',
    'V',
    'W',
    #'X',
    'Y'
)



pairs = []
for i, am in enumerate(amin):
    for am2 in range(i, len(amin)):
        pairs.append(f'{am}{amin[am2]}')




# Lectura de secuencias (placeholder)

seqs = (
    'DADA',
    'AAAC',
    'AACC',
    'AADA',
    'AACC',
    'AADC'
)

# Lectura de secuencias real

filename = input("FASTA a usar: ")


try:
    with open(filename) as file:
    #seqs = [line.rstrip() for line in file if line.rstrip().startswith('>') == False]
    #seqs = [line.rstrip() for line in file]
        seqs = []
        curline = ""
        for line in file:
            #print(line.rstrip())
            if line.rstrip().startswith('>'):
                seqs.append(curline)
                curline = ""
            else:
                curline += line.rstrip()
        seqs.pop(0)
except Exception:
    print('Usando placeholder')

title = input('Title of graph: ')

#print(seqs)


# Paso 1 - contar la frecuencia de cada aminoacido

num_of_amin = np.zeros(len(amin))

length = len(seqs[0])
n_of_seq = len(seqs)

for i, seq in enumerate(seqs):
    for j, aa in enumerate(amin):
        #print(seq.count(aa))
        num_of_amin[j] = num_of_amin[j] + seq.count(aa)
    print(f'\rPaso 1/5: {i} / {n_of_seq}')
total_amin = sum(num_of_amin)

aminDF = pd.DataFrame({
    "Amin": amin,
    "Num": num_of_amin
})

print("Conteo de aminoácidos finalizado")

# Paso 2 - Contar la frecuencia de cada par

num_of_pair = np.zeros(len(pairs))

for coln in range(length):
    curcol = get_col(coln, seqs)
    #print(curcol)
    for j, pair in enumerate(pairs):
        if pair[0] not in curcol or pair[1] not in curcol:
            continue
        #print(pair)
        #N = ( curcol.count(pair[0]) + curcol.count(pair[1]) ) if pair[0] != pair[1] else curcol.count(pair[0])
        if pair[0] == pair[1]:
            N = curcol.count(pair[0])
            num_of_pair[j] += sp.special.comb(N, 2)
        else:
            N_full = ( curcol.count(pair[0]) + curcol.count(pair[1]) )
            N_left = curcol.count(pair[0])
            N_right = curcol.count(pair[1])
            num_of_pair[j] += sp.special.comb(N_full, 2) - sp.special.comb(N_left, 2) - sp.special.comb(N_right, 2)
        #print(num_of_pair[j])
    print(f'\rPaso 2/5: {coln} / {length}')

total_pairs = sum(num_of_pair)
#print(pairs)
#print(num_of_pair)

pairsDF = pd.DataFrame(
    {'Par': pairs,
    'Num': num_of_pair}
)

#print(rs[rs['Par'] == 'AA']['Num'])

#print(rs.to_string())

print("Conteo de pares finalizado")

# Paso 3 - Contar la frecuencia observada de cada par de aminoacidos

freq_pair_obs = num_of_pair / total_pairs
pairsDF["FreqObs"] = freq_pair_obs

# Paso 4 - Contar la frecuencia esperada de cada par

freq_pair_exp = np.zeros(len(freq_pair_obs))

for i, pair in enumerate(pairs):
    #print(float(aminDF[aminDF['Amin'] == pair[0]]['Num']))
    factor = 2 if pair[0] != pair[1] else 1
    freq_pair_exp[i] += ((
        (float(aminDF[aminDF['Amin'] == pair[0]]['Num']) / total_amin)
    * 
        (float(aminDF[aminDF['Amin'] == pair[1]]['Num']) / total_amin)
    ) * factor)
    print(f'\rPaso 4/5: {i} / 230')

pairsDF['FreqExp'] = freq_pair_exp
#print(pairsDF.to_string())

# Paso 5 - Calcular Log-odd ratio

pairsDF['LOR'] = 2 * np.log2(pairsDF['FreqObs'] / pairsDF['FreqExp'])


print("Cálculo de Log odd ratio finalizado")
#print(pairsDF.to_string())

# Paso 6 - Graficas LOR
# X: aa 1
# Y: aa 2
# Z: lor

dtg = pd.DataFrame({
    '1aa': np.concatenate(([pair[0] for pair in pairs], [pair[1] for pair in pairs]), axis=None),
    '2aa': np.concatenate(([pair[1] for pair in pairs], [pair[0] for pair in pairs]), axis=None),
    'LOR': np.concatenate((pairsDF['LOR'], pairsDF['LOR']), axis=None)
})

dtg = dtg.drop_duplicates()

dtgp = dtg.pivot(index='1aa', columns='2aa', values='LOR')

dtgp = dtgp[[
    '-',
    'G',
    'A',
    'V',
    'L',
    'I',
    'S',
    'T',
    'Y',
    'C',
    'M',
    'D',
    'E',
    'N',
    'Q',
    'R',
    'K',
    'H',
    'F',
    'W',
    'P'
]]

dtgp = dtgp.loc[[
    '-',
    'G',
    'A',
    'V',
    'L',
    'I',
    'S',
    'T',
    'Y',
    'C',
    'M',
    'D',
    'E',
    'N',
    'Q',
    'R',
    'K',
    'H',
    'F',
    'W',
    'P'
]]

print(dtg.to_string())

#sns.set_theme()
plt.figure(figsize=(10,6))
ax = sns.heatmap(
    data=dtgp,
    annot=True,
    cmap='inferno',
    linewidth=0.0,
    linecolor='black',
    vmin=-20,
    vmax=20
    )

plt.yticks(rotation=0)

#ax.invert_xaxis()
#ax.invert_yaxis()

plt.xlabel("")
plt.ylabel("")
plt.title(title)
plt.show()


# Legacy

"""
No = 0

A = 0
C = 0
D = 0
E = 0
F = 0
G = 0
H = 0
I = 0
K = 0
L = 0
M = 0
N = 0
P = 0
Q = 0
R = 0
S = 0
T = 0
V = 0
W = 0
X = 0
Y = 0
"""

"""
for i, am in enumerate(amin):
    for am2 in range(i, len(amin)):
        print(f'    "{am}{amin[am2]}",')
"""

"""
for coln in range(length):
    col = ""
    for seq in seqs:
        col + seq[coln]
"""

"""
pairs = (
    "--",
    "-A",
    "-C",
    "-D",
    "-E",
    "-F",
    "-G",
    "-H",
    "-I",
    "-K",
    "-L",
    "-M",
    "-N",
    "-P",
    "-Q",
    "-R",
    "-S",
    "-T",
    "-V",
    "-W",
    "-Y",
    "AA",
    "AC",
    "AD",
    "AE",
    "AF",
    "AG",
    "AH",
    "AI",
    "AK",
    "AL",
    "AM",
    "AN",
    "AP",
    "AQ",
    "AR",
    "AS",
    "AT",
    "AV",
    "AW",
    "AY",
    "CC",
    "CD",
    "CE",
    "CF",
    "CG",
    "CH",
    "CI",
    "CK",
    "CL",
    "CM",
    "CN",
    "CP",
    "CQ",
    "CR",
    "CS",
    "CT",
    "CV",
    "CW",
    "CY",
    "DD",
    "DE",
    "DF",
    "DG",
    "DH",
    "DI",
    "DK",
    "DL",
    "DM",
    "DN",
    "DP",
    "DQ",
    "DR",
    "DS",
    "DT",
    "DV",
    "DW",
    "DY",
    "EE",
    "EF",
    "EG",
    "EH",
    "EI",
    "EK",
    "EL",
    "EM",
    "EN",
    "EP",
    "EQ",
    "ER",
    "ES",
    "ET",
    "EV",
    "EW",
    "EY",
    "FF",
    "FG",
    "FH",
    "FI",
    "FK",
    "FL",
    "FM",
    "FN",
    "FP",
    "FQ",
    "FR",
    "FS",
    "FT",
    "FV",
    "FW",
    "FY",
    "GG",
    "GH",
    "GI",
    "GK",
    "GL",
    "GM",
    "GN",
    "GP",
    "GQ",
    "GR",
    "GS",
    "GT",
    "GV",
    "GW",
    "GY",
    "HH",
    "HI",
    "HK",
    "HL",
    "HM",
    "HN",
    "HP",
    "HQ",
    "HR",
    "HS",
    "HT",
    "HV",
    "HW",
    "HY",
    "II",
    "IK",
    "IL",
    "IM",
    "IN",
    "IP",
    "IQ",
    "IR",
    "IS",
    "IT",
    "IV",
    "IW",
    "IY",
    "KK",
    "KL",
    "KM",
    "KN",
    "KP",
    "KQ",
    "KR",
    "KS",
    "KT",
    "KV",
    "KW",
    "KY",
    "LL",
    "LM",
    "LN",
    "LP",
    "LQ",
    "LR",
    "LS",
    "LT",
    "LV",
    "LW",
    "LY",
    "MM",
    "MN",
    "MP",
    "MQ",
    "MR",
    "MS",
    "MT",
    "MV",
    "MW",
    "MY",
    "NN",
    "NP",
    "NQ",
    "NR",
    "NS",
    "NT",
    "NV",
    "NW",
    "NY",
    "PP",
    "PQ",
    "PR",
    "PS",
    "PT",
    "PV",
    "PW",
    "PY",
    "QQ",
    "QR",
    "QS",
    "QT",
    "QV",
    "QW",
    "QY",
    "RR",
    "RS",
    "RT",
    "RV",
    "RW",
    "RY",
    "SS",
    "ST",
    "SV",
    "SW",
    "SY",
    "TT",
    "TV",
    "TW",
    "TY",
    "VV",
    "VW",
    "VY",
    "WW",
    "WY",
    "YY"
) 
"""

"""
amin = (
    '-',
    'G',
    'A',
    'V',
    'L',
    'I',
    'S',
    'T',
    'Y',
    'C',
    'M',
    'D',
    'E',
    'N',
    'Q',
    'R',
    'K',
    'H',
    'F',
    'W',
    #'X',
    'P'
)
"""

"""
dtgp = dtgp[[
    '-',
    'G',
    'A',
    'V',
    'L',
    'I',
    'S',
    'T',
    'Y',
    'C',
    'M',
    'D',
    'E',
    'N',
    'Q',
    'R',
    'K',
    'H',
    'F',
    'W',
    'P'
]]

dtgp = dtgp.loc[[
    '-',
    'G',
    'A',
    'V',
    'L',
    'I',
    'S',
    'T',
    'Y',
    'C',
    'M',
    'D',
    'E',
    'N',
    'Q',
    'R',
    'K',
    'H',
    'F',
    'W',
    'P'
]]
"""
