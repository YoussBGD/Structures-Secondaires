"""
Ce programme assigne des structures secondaires à des protéines en se basant
sur les coordonnées des atomes dans des fichiers pdb avec la méthode DSSP.
"""


import math
import numpy as np
import pandas as pd
import urllib.request
import os

# ============================================================================
# argparse permet de definir les options pour un programme = ex : -h : aide
# ambertools reduce pour ajouter les hydrogenes
# lineprofiller


class DSSP_assignment:
    def __init__(self, PDB_ID, download, addH, model, chain):
        if download is True:
            self.download_pdb(PDB_ID)

# ============================================================================
    """
    Option :
    Cette fonction permet de récuperer et télecharger le .pdb directement
    depuis le site.
    Le fichier est stocké dans le dosier PDB_file
    """
    def download_pdb(self, file_name):
        path_dl_pdb = 'http://files.rcsb.org/download/{}.pdb'.format(file_name)
        path_output = '../PDB_file/{}.pdb'.format(file_name)
        urllib.request.urlretrieve(path_dl_pdb, path_output)


# ============================================================================
    """
    Option :
    Cette fonction appel ambertools reduce et protone la protéine
    Obligatoire pour les fichiers Xray ne contenant pas les H
    Remplace le fichier .pdb de base
    """
    def add_H(self, file_name):
        os.system(
            "reduce -build -nuclear > /dev/null 2>&1 ../PDB_file/{}.pdb"
            "> ../PDB_file/{}H.pdb".format(file_name, file_name))

        os.system("mv ../PDB_file/{}H.pdb ../PDB_file/{}.pdb"
                  .format(file_name, file_name))
# ============================================================================

    """
    Lit le fichier pdb et stocke ce qui est nécessaire dans un
    dataframe pandas.
    Lit un seul modèle (à choisir) et une seule chaine (à choisir)
    Pour chaque résidus, on retrouve 4 atomes nécessaires pour le calcul DSSP:
    N, H, C, O du backbone
    """
    def read_pdb(self, file_name, model, chain):

        # Récupération du pdb
        file_name = '../PDB_file/{}.pdb'.format(file_name)

        dico_infos_prec = {}
        dico_infos_prec["atom"] = "NULL"

        # Lecture
        with open(file_name, "r") as f_pdb:
            coor_lst = []
            flag = True  # flag pour le choix du modèle

            # Pour renommer les aa à partir de 1
            j = 0
            i = 1

            # Lecture des lignes du pdb
            for ligne in f_pdb:
                # Tous les 4 atomes, on incrémente i
                if j == 4:
                    i += 1
                    j = 0

                # On change le flag en fonction du modèle
                if (ligne.startswith("MODEL")):
                    model_n = ligne.split()
                    if int(model_n[1]) == model:
                        flag = True
                    else:
                        flag = False

                # Confitions pour remplir le dictionnaire
                if (flag is True and ligne.startswith("ATOM") and
                    (ligne[21:23].strip()) in (chain, "") and
                    ((ligne[12:16].strip() in ["H", "H1", "O", "N", "C"]) or
                    (str(ligne[17:21].strip()) == 'PRO' and
                        str(ligne[12:16].strip()) == "HA"))):

                    # Création du dictionnaire vide.
                    dico_infos = {}

                    # Extraction de l'atome
                    dico_infos["atom"] = str(ligne[12:16].strip())

                    # Extraction du nom du résidu.
                    dico_infos["residu"] = str(ligne[17:21])

                    # aa renommés à partir de 1 pour éviter certains problèmes
                    dico_infos["N° resid"] = i

                    # Extraction de la coordonnée x.
                    dico_infos["X"] = float(ligne[30:38])

                    # Extraction de la coordonnée y.
                    dico_infos["Y"] = float(ligne[38:46])

                    # Extraction de la coordonnée z.
                    dico_infos["Z"] = float(ligne[46:54])

                    # Extraction de la chaine
                    dico_infos["chaine"] = str(ligne[21:23])

                    # Condition qui permet d'éviter certains problèmes
                    # dûs à ambertools_reduce qui peut ajouter plusieurs fois
                    # le même H = on ne garde que le premier
                    # pour toujours avoir 4 atomes / residus
                    if (dico_infos["atom"] != dico_infos_prec["atom"]):
                        coor_lst.append(dico_infos)
                        j = j+1

                    # Pour pouvoir comparer si 2H d'affilé
                    dico_infos_prec["atom"] = dico_infos["atom"]

            # DF rempli à chaque ittération
            self.data_frame = pd.DataFrame(coor_lst)

# ============================================================================

    """
    Permet de retrouver un résidu d'intérêt rapidement
    Nécessite qu'on ne possède que 4 atomes / résidus et que les
    numéros des résidus commence à 1 et se suivent
    """
    def coord_un_res(self, num_res):
        dic_coord = {}
        dico_atoms = {}
        liste_atom = []
        if num_res == 1:
            i = 0
        else:
            # pour accéder directement au res d'intérêt sans tout parcourir
            # on fait d'un atome on fait (numero du résidu-1)*4
            i = 4 * (num_res - 1)

        while self.data_frame["N° resid"][i] == num_res:

            dic_coord = {  # Dictionnaire de coordonnees
                    "X": float(self.data_frame["X"][i]),
                    "Y": float(self.data_frame["Y"][i]),
                    "Z": float(self.data_frame["Z"][i])}

            dico_atoms[str(self.data_frame["atom"][i])] = dic_coord
            dic_coord = {}
            i = i + 1
        return dico_atoms


# ============================================================================

    """
    Calcul des ditances
    """
    def dist_3D(self, A, B):
        dst3d = math.sqrt((A["X"] - B["X"]) ** 2 + (A["Y"] - B["Y"]) ** 2 +
                          (A["Z"] - B["Z"]) ** 2)
        return dst3d

# ============================================================================

    """
    Calcul des énergies entre 2 résidus selon l'algorithme DSSP
    Fait appel à dist_3D
    """
    def calc_Energie(self, A, B):
        Q1_CO = 0.42
        Q2_NH = 0.20

        # Cas normal
        if("H" in B.keys()):
            E = Q1_CO * Q2_NH * ((1 / self.dist_3D(B["N"], A["O"])) +
                                 (1 / self.dist_3D(B["H"], A["C"])) -
                                 (1 / self.dist_3D(B["H"], A["O"])) -
                                 (1 / self.dist_3D(B["N"], A["C"]))) * 332

        # Autres cas :
        # N-ter : on a choisit H1 pour remplir le Df (pour avoir 4 atomes/res)
        # Proline : pas de H sur le N --> on a choisit HA pour remplir le DF
        # Pour ces 2 cas, on met E=0 car impossibilité de liaisons H DSSP

        else:
            E = 0

        return E

# ============================================================================

    """
    Calcul de la matrice d'energie
    Fait appel à coord_un_res pour trouver les 2 residus
    Puis appel calc_Energie pour calculer l'énergie entre 2 résidus
    """
    def calc_matrice_energie(self):
        length = int(len(self.data_frame) / 4)
        matrice = np.zeros((length, length))
        for i in range(1, length):
            A = self.coord_un_res(i)
            for j in range(i+1, length):
                B = self.coord_un_res(j)
                E = self.calc_Energie(A, B)
                matrice[i-1, j-1] = E
                matrice[j-1, i-1] = E

        self.matrice_energie = matrice

# ============================================================================
    """
    Assigne les hélices alpha selon DSSP
    """
    def helice_alpha(self, ss_elements, type):
        alpha = []

        # Helice 3,4 ou 5
        for i in range(0, (self.matrice_energie.shape[0] - type)):
            j = i + type
            if (ss_elements[i] == "C" and ss_elements[j] == "C"):
                E = self.matrice_energie[i, j]
                if E < - 0.5:
                    alpha.append([i + 1, j + 1])
            # alpha = liste de liste contenant les résidus i-j en interaction

            # Transformation en array
            alpha_array = np.asarray(alpha)

            # Permet de passer de 2D à 1D
            # Ex : alpha_array= [[i, i+4] ,[i+1 , i+5]] --> [i, i+4, i+1, i+5]
            res_alpha = np.ravel(alpha_array)

            # Garde que les uniques
            # On pouvait avoir des doublons :
            # Par exemple [i; i+4] et [i+4; i+8] donc i+4 présent 2 fois
            res_alpha = np.unique(res_alpha)
            res_alpha = np.sort(res_alpha)
        for i in res_alpha:
            if type == 4:
                ss_elements[i] = "H"
            if type == 3:
                ss_elements[i] = "G"
            if type == 5:
                ss_elements[i] = "I"

        # Regulation des hélices puis return
        ss_elements = self.regulation_ss(ss_elements)
        return ss_elements

# ============================================================================
    """
    Assigne les feuillets B selon DSSP
    44 aa minimum pour la séparation de 2 feuillets car on a au moins un coude
    de 3 aa
    """
    def feuillet_beta(self, ss_elements):

        pont_PB = []
        pont_APB = []
        list_PB = []
        list_APB = []
        for i in range(1, (self.matrice_energie.shape[0] - 1)):
            for j in range(i + 1, (self.matrice_energie.shape[0] - 1)):
                if (ss_elements[i] == "C" and ss_elements[j] == "C" and
                        abs(i - j) > 4):

                    # Parallele
                    E_PB1 = [self.matrice_energie[i - 1, j],
                             self.matrice_energie[j, i + 1]]

                    E_PB2 = [self.matrice_energie[j - 1, i],
                             self.matrice_energie[i, j + 1]]

                    if ((E_PB1[0] < -0.5 and E_PB1[1] < - 0.5) or
                            (E_PB2[0] < -0.5 and E_PB2[1] < - 0.5)):

                        pont_PB = [i + 1, j + 1]
                        list_PB.append(pont_PB)

                    pont_PB = []

                    # Antiparallele
                    E_APB1 = [self.matrice_energie[i, j],
                              self.matrice_energie[j, i]]

                    E_APB2 = [self.matrice_energie[i - 1, j + 1],
                              self.matrice_energie[j - 1, i + 1]]

                    if ((E_APB1[0] < - 0.5 and E_APB1[1] < - 0.5) or
                            (E_APB2[0] < - 0.5 and E_APB2[1] < - 0.5)):
                        pont_APB = [i + 1, j + 1]
                        list_APB.append(pont_APB)

                    pont_APB = []

        PB_array = np.asarray(list_PB)
        res_PB = np.ravel(PB_array)
        res_PB = np.unique(res_PB)
        res_PB = np.sort(res_PB)
        res_PB = res_PB

        for i in res_PB:
            ss_elements[i] = "E"

        APB_array = np.asarray(list_APB)
        res_APB = np.ravel(APB_array)
        res_APB = np.unique(res_APB)
        res_APB = np.sort(res_APB)
        res_APB = res_APB

        for i in res_APB:
            ss_elements[i] = "E"

        # Regulation des feuillets puis return
        ss_elements = self.regulation_ss(ss_elements)
        return ss_elements

# ============================================================================
    def beta_turn(self, ss_elements, type):
        turn = []
        res_turn = []
        # turn 3,4 ou 5 / les turns sont un ensemble de residus dont deux
        # sont reliés avec une liaison H
        for i in range(0, (self.matrice_energie.shape[0]-type)):
            j = i + type

            # test si c'est un coil et si i ou j n'ont pas été deja mis dans
            # la liste pour ne pas avoir des doublons.
            if (ss_elements[i] == "C" and ss_elements[j] == "C" and
                    i not in res_turn and j not in res_turn):

                E = self.matrice_energie[i, j]
                if E < - 0.5:
                    turn.append(list(range(i, j)))
                    turn_array = np.asarray(turn)
                    res_turn = np.ravel(turn_array)
        res_turn = np.unique(res_turn)
        res_turn = np.sort(res_turn)

        for i in res_turn:
            ss_elements[i] = "T"

        return ss_elements

# ============================================================================

    """
    La fonction suivante verifie si les structures secondaires respectent bien
    leurs regles , par exemple elle teste si un brin beta ou une helice alpha
    n est pas constitue d un seul acide amine etc..
    """

    def regulation_ss(self, ss_elements):
        ss_elements2 = ''.join(ss_elements)
        ss = []
        deb = 0
        fin = 0
        for i in range(len(ss_elements2)):
            if(i != len(ss_elements2) - 1 and
               ss_elements2[i] != ss_elements2[i + 1]):

                fin = i
                liste = ss_elements2[deb:fin+1]
                deb = fin + 1
                ss.append(liste)
        ss.append("".join(ss_elements2[deb:len(ss_elements2)]))
        ss2 = list(ss)
        for c in range(len(ss)):
            if(len(ss2[c]) == 1 and c != len(ss) - 1):
                ss2[c] = ss2[c + 1][0]
            if(len(ss2[c]) == 1 and c == len(ss) - 1):
                ss2[c] = ss2[c - 1][0]
        for struct in range(len(ss2)):
            if(ss2[struct][0] == "G" and len(ss2[struct]) < 3):
                ss2[struct] = "C" * len(ss2[struct])
            if(ss2[struct][0] == "H" and len(ss2[struct]) < 4):
                ss2[struct] = "C" * len(ss2[struct])
            if(ss2[struct][0] == "I" and len(ss2[struct]) < 5):
                ss2[struct] = "C" * len(ss2[struct])
            if(ss2[struct][0] == "T" and len(ss2[struct]) < 3):
                ss2[struct] = "C" * len(ss2[struct])

        ss_elements = "".join(ss2)
        ss_elements = list(ss_elements)

        return ss_elements
# ============================================================================
    """
    Calcul des SS :
    - De base : toute la protéine est assignée en random coil : liste de "C"
    - Ensuite, les hélices sont assignées (modification des "C" en "H")
    - Puis les feuillets sont assignées la ou il n'y a pas de H ("C" en "E")
    et ainsi de suite
    - Ordre d'assignation DSSP : H > E > G > I > T
    - A chaque étape on régule pour éviter d'avoir des structures secondaires
    d'un seul acide aminé, ou un seul cassant une hélice...
    """
    def calc_ss(self):
        ss_elements = list("C"*(self.matrice_energie.shape[0]))
        ss_elements = self.helice_alpha(ss_elements, 4)  # Calc_ss+ regulation
        ss_elements = self.feuillet_beta(ss_elements)  # Calc_ss + regulation
        ss_elements = self.helice_alpha(ss_elements, 3)  # Calc_ss+ regulation
        ss_elements = self.helice_alpha(ss_elements, 5)  # Calc_ss+ regulation
        # Assignation de tous les turns puis regulation
        ss_elements = self.beta_turn(ss_elements, 3)
        ss_elements = self.beta_turn(ss_elements, 4)
        ss_elements = self.beta_turn(ss_elements, 5)
        ss_elements = self.regulation_ss(ss_elements)
        self.ss_elements = ss_elements

# ============================================================================

    # Calculs de pourcentage des structures secondaires dans une seule
    # frame ou dans une trajectoire entiere et ecriture dans un fichier
    def pourcentage_ss(self):
        dic_ss = {}
        dic_ss["C"] = (self.ss_elements.count("C") /
                       len(self.ss_elements)) * 100
        dic_ss["G"] = (self.ss_elements.count("G") /
                       len(self.ss_elements)) * 100
        dic_ss["I"] = (self.ss_elements.count("I") /
                       len(self.ss_elements)) * 100
        dic_ss["H"] = (self.ss_elements.count("H") /
                       len(self.ss_elements)) * 100
        dic_ss["T"] = (self.ss_elements.count("T") /
                       len(self.ss_elements)) * 100
        dic_ss["E"] = (self.ss_elements.count("E") /
                       len(self.ss_elements)) * 100
        return dic_ss

    def write_pourcentage(self, list_dic, PDB_ID):
        with open("../Output/Stat_{}.txt".format(PDB_ID), "w") as Output:
            i = 1
            for frame in list_dic:
                Output.write("Frame  {}\n".format((i)))
                for ss in frame:
                    if ss == "C":
                        Output.write("Cette protéine contient {:.2f} %  de "
                                     "Coils\n".format((frame[ss])))
                    if ss == "G":
                        Output.write("Cette protéine contient {:.2f} %  de "
                                     "Right-handed 3_10 helix\n"
                                     .format((frame[ss])))
                    if ss == "H":
                        Output.write("Cette protéine contient {:.2f} %  de "
                                     "Right-handed α-helix\n"
                                     .format((frame[ss])))
                    if ss == "I":
                        Output.write("Cette protéine contient {:.2f} %  de "
                                     "Right-handed π-helix\n"
                                     .format((frame[ss])))
                    if ss == "T":
                        Output.write("Cette protéine contient {:.2f} %  de "
                                     "Turns\n"
                                     .format((frame[ss])))
                    if ss == "E":
                        Output.write("Cette protéine contient {:.2f} %  de "
                                     "strand β-pleated sheet\n"
                                     .format((frame[ss])))

                Output.write("\n\n")
                i = i + 1
