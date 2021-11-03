"""
Ce script permet d'afficher les resultats d'assignations des structures
secondaires sur une trajectoire -> Residus en fonction des frames
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def DSSP_trj_plot(liste):

	# Changement des lettres du dictionnaire DSSP pour les legendes
    for frame in liste:

        for res in range(len(frame)):

            if frame[res] == "H":
                frame[res] = "A-helix"
            if frame[res] == "E":
                frame[res] = "B-Sheet"
            if frame[res] == "G":
                frame[res] = "3-Helix"
            if frame[res] == "I":
                frame[res] = "Pi-Helix"
            if frame[res] == "C":
                frame[res] = "Coil"
            if frame[res] == "T":
                frame[res] = "Turn"

    df = pd.DataFrame(liste)
    df = df.transpose()
    # Pour commencer au résidus 1 et frame 1 et pas à 0
    df.columns = list(range(1, len(liste) + 1))
    df.index = list(range(1, len(liste[0]) + 1))

    value_to_int = {j: i for i, j in
                    enumerate(pd.unique(df.values.ravel()))}

    n = len(value_to_int)

    cmap = sns.color_palette("bright", n)

    Xframe_axis_labels = list(range(1, len(liste) + 1))
    Yres_axis_labels = list(range(1, len(liste[0]) + 1))

    ax = sns.heatmap(df.replace(value_to_int), cmap=cmap,
                     xticklabels=5, yticklabels=10)

    colorbar = ax.collections[0].colorbar
    r = colorbar.vmax - colorbar.vmin
    colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
    colorbar.set_ticklabels(list(value_to_int.keys()))
    ax.invert_yaxis()
    ax.set_xlabel("Frame")
    ax.set_ylabel("Residus")

    plt.show()
