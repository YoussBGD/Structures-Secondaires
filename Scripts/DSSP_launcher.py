"""
Ce programme appel les scripts Draw_ss.py, Draw_ss_trajectory.py et
DSSP_functions.py et gère les différentes options décrites ci-dessous.
Il s'agit du script à executer pour utiliser notre programme DSSP
"""

import argparse
import DSSP_functions
import Draw_ss_trajectory
import Draw_ss
import pandas as pd
import numpy as np
import sys
import time


# ============================================================================
"""
Toutes les options proposées à l'utilisateur
"""
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("pdb_file", help="Give PDB ID - Exemple : 4o8h")
parser.add_argument("-m", "--model", type=int, default=1,
                    help="Model number (for RMN structure) - Default = 1")
parser.add_argument("-c", "--chain", default="A",
                    help="Chain number to analyze - Default = A")
parser.add_argument("-d", "--downloadPDB", action="store_true",
                    help="Download file from PDB - Default = False")
parser.add_argument("-nh", "--no_addH", action="store_true",
                    help="Don't use Ambertools to add H - Default = False")
parser.add_argument("-p", "--plotSS", action="store_true",
                    help="Secondary structure plot - Default = False")
parser.add_argument("-trj", "--DSSP_Trajectory", action="store_true",
                    help="Do DSSP on a trajectory (on PDB file)"
                    " - Default = False")
parser.add_argument("-stat", "--statistiques", action="store_true",
                    help="To calculate statistics of Secondary structures"
                    "(on PDB file) - Default = False")


"""
Récupération des options
"""
args = parser.parse_args()
pdb_file = args.pdb_file
model = args.model
chain = args.chain
no_add_H = args.no_addH
plot_SS = args.plotSS
download_PDB = args.downloadPDB
traj = args.DSSP_Trajectory
stat = args.statistiques

# ============================================================================

"""
Barre de progression pour les fichiers de trajectroires
"""


def progressbar(it, prefix="", size=60, file=sys.stdout):

    count = len(it)

    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] Frame %i/%i\r" %
                   (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()

    show(0)

    for i, item in enumerate(it):
        yield item
        show(i+1)

    file.write("\n")
    file.flush()

# ============================================================================

if __name__ == "__main__":

    """
    1 fichier pdb
    """
    if traj is False:
        # Appel de la classe du script DSSP_functions.py
        monpdb = DSSP_functions.DSSP_assignment(pdb_file, download_PDB,
                                                no_add_H, model, chain)

        # Appel à la fonction add_H de la classe
        if no_add_H is False:
            print("\nProtonation par Ambertools reduce en cours...")
            monpdb.add_H(pdb_file)
            print("Protonation terminée\n")

        # Appel à la fonction read_pdb de la classe
        monpdb.read_pdb(pdb_file, model, chain)

        # Appel à la fonction calculant la matrice d'énergie
        # entre les résidus de la classe
        print("Calcul de la matrice d'énergie")
        monpdb.calc_matrice_energie()

        # Appel à la fonctiondf=pd.DataFrame(liste) calculant
        # les structures secondaires de la classe
        print("Assignation des structures secondaires \n")
        monpdb.calc_ss()

        # Récupération de la liste contenant la SS calculé pour chaque résidus
        ss_elements = monpdb.ss_elements

        # Ecriture de l'assignation dans un fichier de sortie
        with open("../Output/SS_{}.txt".format(pdb_file), "w") as f:
            for i in range(len(ss_elements)):
                f.write("Résidu {} : {}\n".format((i+1), ss_elements[i]))

        # Options de sorties

        # Ecriture des statistiques
        prc_ss_element = [monpdb.pourcentage_ss()]
        if stat is True:
            monpdb.write_pourcentage(prc_ss_element, pdb_file)

        # Plot des SS (a helix, b strand) si voulu
        if plot_SS is True:
            Draw_ss.visualize_secondary_structure(ss_elements, 1)

    """
    Trajectory
    """
    if traj is True:
        # Récupération du nombre de frame
        all_frame = 0
        with open("../PDB_file/{}.pdb".format(pdb_file), "r") as file:
            for line in file:
                if (line.startswith("MODEL")):
                    model = line.split()
                    all_frame = int(model[1])

        # Appel de la classe du script DSSP_functions.py
        monpdb = DSSP_functions.DSSP_assignment(pdb_file, download_PDB,
                                                no_add_H, model, chain)

        # Appel à la fonction add_H de la classe
        if no_add_H is False:
            print("\nProtonation par Ambertools reduce en cours...")
            monpdb.add_H(pdb_file)
            print("Protonation términée\n")

        # Boucle sur toutes les frames (MODEL 1 --> MODEL X)
        # Stockage des résultats dans une liste de liste
        trj_ss_elements = []
        trj_prc_ss_element = []
        print("\nDébut des calculs :\n")
        for frame in progressbar(range(1, all_frame+1), "Computing: ", 50):
            monpdb.read_pdb(pdb_file, frame, chain)
            monpdb.calc_matrice_energie()
            monpdb.calc_ss()
            ss_elements = monpdb.ss_elements
            trj_ss_elements.append(ss_elements)
            prc_ss_element = monpdb.pourcentage_ss()
            trj_prc_ss_element.append(prc_ss_element)
            # Statistiques en sortie pour chaque frame
            if stat is True:
                monpdb.write_pourcentage(trj_prc_ss_element, pdb_file)

        # Ecriture des assignations de chaque frame
        # dans un fichier de sortie csv
        df = pd.DataFrame(trj_ss_elements)
        df.columns = ["Res " + str(i) for i in range(1,
                      len(ss_elements)+1)]
        df.index = ["Frame " + str(i) for i in range(1,
                    len(trj_ss_elements)+1)]
        output = "../Output/SS_{}.csv".format(pdb_file)
        df.to_csv(output, header=True)

        # Option de sortie
        # Plot des SS au cours du temps sous forme de heatmap
        if plot_SS is True:
            Draw_ss_trajectory.DSSP_trj_plot(trj_ss_elements)
