# Sorbonne Université 3I024 2023-2024
# TME 2 : Cryptanalyse du chiffre de Vigenere
#
# Etudiant.e 1 : SRAIEB 21114754
# Etudiant.e 2 : SOUAIBY 21102782

import sys, getopt, string, math , collections
from decimal import Decimal

# Alphabet français
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

# Fréquence moyenne des lettres en français
# À modifier
freq_FR = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#resultat de la modification de freq_FR, tableau de reference a la langue francaise : freq_FR = freq(read("germinal.txt"))
freq_FR = [72501.0, 8148.0, 23748.0, 29538.0, 135149.0, 8608.0, 8353.0, 8434.0, 59075.0, 3016.0, 55.0, 48291.0, 20852.0, 55326.0, 38669.0, 18648.0, 7995.0, 52009.0, 61511.0, 58029.0, 50017.0, 12945.0, 9.0, 3204.0, 1810.0, 965.0]

# Chiffrement César
def chiffre_cesar(txt, key):
    """
    Chiffre un texte en utilisant le chiffrement de César.
    Prend en argument le texte à chiffrer et la clé de chiffrement (un entier représentant le décalage dans l'alphabet).
    Retourne texte chiffré (décalé).
    """

    txt_chiffre = ""
    for char in txt:
        if char.isalpha(): #verifier que c est une lettre
            shift = 65 if char.isupper() else 97 #code ascii des caracteres respectifs a et A
            txt_chiffre += chr((ord(char) - shift + key) % 26 + shift) #retrouve le code ascii apres decalage puis on reconvertit en caractere
        else:
            txt_chiffre += char
    return txt_chiffre

# Déchiffrement César
def dechiffre_cesar(txt, key):
    """
    Déchiffre un texte chiffré avec le chiffrement de César.
    """
    return chiffre_cesar(txt, -key) #decalage inverse pour retrouver le texte initial

# Chiffrement Vigenere

def chiffre_vigenere(txt, key):
    """
    Chiffre un texte avec la méthode de chiffrement de Vigenère.
    Prend en argument le texte à chiffrer et la clé de chiffrement (liste d' entiers représentant des décalages dans l'alphabet).
    Retourne texte chiffré.
    """
    txt_chiffre = ''
    n = len(key)
    idx_key = 0
    for char in txt:
        if char.isalpha(): #verifier que c est une lettre 
            txt_chiffre += chiffre_cesar(char, key[idx_key])  # Chiffrement de chaque caractère avec le décalage correspondant
            idx_key = (idx_key + 1) % n  # Passage à la lettre suivante de la clé
        else:
            txt_chiffre += char  # Conserver les caractères qui ne sont pas des lettres
    return txt_chiffre


# Déchiffrement Vigenere
def dechiffre_vigenere(txt, key):
    """
    Déchiffre un texte chiffré avec la méthode de chiffrement de Vigenère.
    """
    #meme principe que le code precedent sauf pour cesar : dechiffrement au lieu de chiffrement
    txt_dechiffre = ''
    n = len(key)
    idx_key = 0 
    for char in txt:
        if char.isalpha():
            txt_dechiffre += dechiffre_cesar(char, key[idx_key]) 
            idx_key = (idx_key + 1) % n 
        else:
            txt_dechiffre += char 
    return txt_dechiffre


# Analyse de fréquences

def freq(txt):
    """
    Calcule les fréquences d'apparition des lettres dans un texte.
    Prend en argument le texte pour lequel calculer les fréquences.
    Retourne une liste contenant les fréquences d'apparition de chaque lettre dans le texte.
    """
    txt = txt.upper()  # Convertir le texte en majuscule pour normaliser si le texte n'est pas nettoyé
    frequencies = {lettre: 0 for lettre in alphabet}
    frequencies.update(collections.Counter(char for char in txt if char.isalpha())) # dictionnaire contenant les nombres d'occurences des lettres du texte.
    ordered_frequencies = { lettre: frequencies.get(lettre, 0) for lettre in sorted(frequencies) }  # trier par ordre alphabetique
    return [val * 1.0 for val in ordered_frequencies.values()]


# Renvoie l'indice dans l'alphabet
# de la lettre la plus fréquente d'un texte
def lettre_freq_max(txt):
    """
    Renvoie l'indice dans l'alphabet de la lettre la plus fréquente dans un texte.
    Prend en argument le texte pour lequel trouver la lettre la plus fréquente.
    Retourne l'indice dans l'alphabet (de 0 à 25) de la lettre la plus fréquente dans le texte.
    """
    frequences = freq(txt)  # Calcul des fréquences d'apparition des lettres dans le texte
    max = 0
    imax = 0
    for i in range(len(frequences)):
        if frequences[i]>max:
            imax = i
            max = frequences[i]  
            # Recherche de l'indice de la lettre la plus fréquente
    return imax


# indice de coïncidence
def indice_coincidence(hist):
    """
    Calcule l'indice de coïncidence d'un texte.
    Prend en argument un tableau qui correspond aux occurrences des lettres d'un texte.
    retourne l'indice de coïncidence du texte.
    """
    total_chars = sum(hist)  # Nombre total de lettres dans le texte
    ic_sum = sum(ni * (ni - 1) for ni in hist)  # Somme des produits ni * (ni - 1) avec ni le nombre d'occurrences de la lettre d'indice i dans le texte.

    if total_chars <= 1:  # Si le texte est vide ou a seulement une lettre
        return 0.0

    ic = ic_sum / (total_chars * (total_chars - 1))  # Calcul de l'indice de coïncidence selon la formule
    return ic


# Recherche la longueur de la clé
def longueur_clef(cipher):
    """
    Détermine la longueur probable de la clé utilisée pour chiffrer le texte.
    """
    max_key_length = 20 
    avg_ic_threshold = 0.06  # Seuil de l'indice de coincidence moyen pour considérer une taille de clé

    for key_length in range(1, max_key_length + 1):
        
        columns = [cipher[i::key_length] for i in range(key_length)] # Découper le texte en colonnes en utilisant la taille de clé actuelle
        avg_ic = sum(indice_coincidence(freq(column)) for column in columns) / key_length #indice de coïncidence moyen pour chaque colonne

        if avg_ic > avg_ic_threshold:
            return key_length

    return 0

    
# Renvoie le tableau des décalages probables étant
# donné la longueur de la clé
# en utilisant la lettre la plus fréquente
# de chaque colonne
def clef_par_decalages(cipher, key_length):
    """
    Détermine la clé probable utilisée pour chiffrer le texte avec la méthode de Vigenère.
    Prend le texte chiffré et la longueur de la clé, retourne la clé sous forme d'une table de décalages.
    """
    key = []

    for i in range(key_length): # iterer sur chaque colonne
        column = cipher[i::key_length]
        most_frequent_letter_index = lettre_freq_max(column) #indice de la lettre la plus fréquente dans la colonne
        shift = (most_frequent_letter_index - 4) % 26 #décalage par rapport à la lettre 'E' (indice 4 dans l alphabet)
        key.append(shift)
    return key


# Cryptanalyse V1 avec décalages par frequence max
def cryptanalyse_v1(cipher):
    """
    Effectue une première forme de cryptanalyse sur un texte chiffré.
    Prend le texte chiffré et retourne un texte déchiffré par une clé trouvée par la méthode de l'indice de coincidence.
    """
    key_length = longueur_clef(cipher)
    if key_length == 0: #si on n'arrive pas à déterminer la taille de la clé on ne peut pas avancer
        return cipher

    key = clef_par_decalages(cipher, key_length)

    return dechiffre_vigenere(cipher, key)

#18 textes sont correctement cryptanalysés parmi 100, soit moins de 20% des textes. Cela s'explique notamment par le choix de la longueur de la clé.
#En effet notre algorithme choisit la première longueur qui dépasse le seuil de la moyenne d'indice de coincidence, mais sans aucune garantie qu'elle soit la bonne,
#car on pourrait retrouver pour des clés de taille plus grande un meilleur indice.



################################################################


### Les fonctions suivantes sont utiles uniquement
### pour la cryptanalyse V2.

# Indice de coincidence mutuelle avec décalage
def indice_coincidence_mutuelle(h1, h2, d):
  """
  Calcule l'indice de coïncidence mutuelle entre deux textes décalés de d positions.
  Prend en parametre les tableaux des fréquences des lettres dans les 2 textes ainsi que le decalage du second et retourne l'indice de coincidence correspondant.
  """

  n1 = sum(h1) # nbre total de lettres dans chaque texte
  n2 = sum(h2)

  icm = 0
  for i in range(26):
    icm += h1[i] * h2[(i + d) % 26] #en chiffrant le 2e texte avec cesar, la frequence de la lettre cherchee devient celle de la lettre decalee de d positions.

  return icm / (n1 * n2)


# Renvoie le tableau des décalages probables étant
# donné la longueur de la clé
# en comparant l'indice de décalage mutuel par rapport
# à la première colonne


def tableau_decalages_ICM(cipher, key_length):
    """
    Calcule le décalage relatif de chaque colonne par rapport à la première colonne en utilisant l'indice de coïncidence mutuelle.
    Prend le texte chiffré et une longueur de clef et retourne le tableau des decalages optimaux.
    """

    text_length = len(cipher)
    shifts = [0] * key_length

    for i in range(1, key_length):  # Parcours des colonnes
        
        max_icm = 0
        optimal_shift = 0
        for d in range(26):
            icm = indice_coincidence_mutuelle(freq(cipher[::key_length]), freq(cipher[i::key_length]), d) # calcul icm pour chaque décalage
            if icm > max_icm:
                max_icm = icm
                optimal_shift = d

        shifts[i] = optimal_shift
    return shifts





# Cryptanalyse V2 avec décalages par ICM
def cryptanalyse_v2(cipher):
    """
    Effectue une cryptanalyse du texte chiffré en utilisant la méthode de cryptanalyse basée sur l'indice de coincidence mutuelle.
    Prend en parametre le texte chiffré à analyser et retourne un texte clair.
    """
    key_length = longueur_clef(cipher)
    plain_list = list(cipher)  # Convertir la chaine en liste de caractères
    decalages = tableau_decalages_ICM(cipher, key_length)
    decalages[0] = lettre_freq_max(cipher[::key_length]) - 4  # 4 étant l'indice de E, correspondant à la lettre la plus fréquente du chiffré
    for i in range(1, key_length):
        decalages[i] += decalages[0]
    
    for col in range(key_length):
        plain_list[col::key_length] = dechiffre_cesar(cipher[col::key_length], decalages[col]) # dechiffrement César pour chaque colonne

    plain = ''.join(plain_list) #reconversion de la liste en chaine de caracteres
    
    return plain

#38 textes parmi 100 sont correctement cryptanalysés (soit 38%). Cela s'explique par le choix du décalage de la première colonne, plus specifiquement au niveau du choix de la lettre qui correspond à E. 
#En effet, il est fort probable que la lettre la plus fréquente de ce bout de texte et la lettre la plus frequente de la langue (en l occurence E) coincident, mais ce n'est pas necessairement toujours le cas. Il faudrait une analyse un peu plus profonde à ce niveau.


################################################################


### Les fonctions suivantes sont utiles uniquement
### pour la cryptanalyse V3.

# Prend deux listes de même taille et
# calcule la correlation lineaire de Pearson
def correlation(x, y):
  """
  Calcule la corrélation de Pearson entre deux listes de même longueur.
  """
  #utilisation de Decimal par souci de precision des valeurs

  x_mean = Decimal(sum(x)) / Decimal(len(x)) #moyennes des 2 listes
  y_mean = Decimal(sum(y)) / Decimal(len(y))

  x_var = Decimal(sum((Decimal(x_i) - x_mean)**2 for x_i in x)) / Decimal(len(x)) #variances des 2 listes
  y_var = Decimal(sum((Decimal(y_i) - y_mean)**2 for y_i in y)) / Decimal(len(y))

  cov_xy = Decimal(sum((Decimal(x_i) - x_mean) * ( Decimal(y_i) - y_mean) for x_i, y_i in zip(x, y))) / Decimal(len(x)) #covariance des 2 listes

  correlation = cov_xy / (Decimal(x_var).sqrt() * Decimal(y_var).sqrt())
  
  return float(correlation)


# Renvoie la meilleur clé possible par correlation
# étant donné une longueur de clé fixée
def clef_correlations(cipher, key_length):
  """
  Calcule la clé qui maximise la corrélation entre chaque colonne du texte chiffré et un texte français,
  et renvoie la moyenne des corrélations maximales obtenues ainsi que la clé sous forme de tableau de décalages.
  """

  #initialisation tableau des decalages et des correlations optimales
  shifts = [0] * key_length
  correlations = [0] * key_length

  for i in range(key_length):  # Parcours des colonnes
        max_corr = 0
        optimal_shift = 0
        for d in range(26):
            freq_col = freq(cipher[i::key_length])
            freq_shifted_col = [freq_col[(i+d)%26] for i in range(26)] #frequences des lettres decalees de la colonne
            corr = correlation(freq_shifted_col, freq_FR) # calcul la correlation de la colonne avec le texte de reference pour chaque décalage
            if corr > max_corr:
                max_corr = corr
                optimal_shift = d

        shifts[i] = optimal_shift
        correlations[i] = max_corr
  
  moy_correlations = sum(correlations)/len(correlations) if len(correlations) > 0 else 0.0
  return (moy_correlations, shifts)



# Cryptanalyse V3 avec correlations
def cryptanalyse_v3(cipher):
  """
  Effectue une cryptanalyse du texte chiffré en utilisant la méthode de cryptanalyse basée sur la corrélation de Pearson.
  Prend en parametre le texte chiffré à analyser et retourne un texte clair.
  """
  max_moy_corr = 0.0
  opt_shifts = []
  for key_length in range(1, 21):
    moy_corr, shifts = clef_correlations(cipher, key_length)
    if moy_corr > max_moy_corr:
        max_moy_corr = moy_corr
        opt_shifts = shifts

  plain_text = dechiffre_vigenere(cipher, opt_shifts)
  return plain_text

#94 textes sont correctement cryptanalysés parmi 100, soit 94%. 
#les 6 textes qui echouent sont les textes 81, 86, 88, 89, 94 et 96.
#les textes en question se caracterisent par la presence de certains mots et expressions en familier (gars, pagaille, pizza, ...) ou peu communs (creationnisme, astrophysicien, ...) ce qui pourrait affecter leur correlation avec notre texte de reference de Zola ecrit dans un langage plus soutenu. 





################################################################
# NE PAS MODIFIER LES FONCTIONS SUIVANTES
# ELLES SONT UTILES POUR LES TEST D'EVALUATION
################################################################


# Lit un fichier et renvoie la chaine de caracteres
def read(fichier):
    f=open(fichier,"r")
    txt=(f.readlines())[0].rstrip('\n')
    f.close()
    return txt

# Execute la fonction cryptanalyse_vN où N est la version
def cryptanalyse(fichier, version):
    cipher = read(fichier)
    if version == 1:
        return cryptanalyse_v1(cipher)
    elif version == 2:
        return cryptanalyse_v2(cipher)
    elif version == 3:
        return cryptanalyse_v3(cipher)

def usage():
    print ("Usage: python3 cryptanalyse_vigenere.py -v <1,2,3> -f <FichierACryptanalyser>", file=sys.stderr)
    sys.exit(1)

def main(argv):
    size = -1
    version = 0
    fichier = ''
    try:
        opts, args = getopt.getopt(argv,"hv:f:")
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt == '-h':
            usage()
        elif opt in ("-v"):
            version = int(arg)
        elif opt in ("-f"):
            fichier = arg
    if fichier=='':
        usage()
    if not(version==1 or version==2 or version==3):
        usage()

    print("Cryptanalyse version "+str(version)+" du fichier "+fichier+" :")
    print(cryptanalyse(fichier, version))
    
if __name__ == "__main__":
   main(sys.argv[1:])
