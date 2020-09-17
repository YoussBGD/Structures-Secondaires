import math


def lecture_fichier_pdb(file_name):
    with open(file_name,"r") as file:
        num_atom={}
        liste_atoms=[]
        dico_atoms={}
        for ligne in file:
            if ligne.startswith("ATOM") and ligne[12:16].strip() in ["H","H1","H2","H3","O","N","C"]:
                dico_coord={  #Dictionnaire de coordonnees
                "x": float(ligne[30:38].strip()), 
                "y": float(ligne[38:46].strip()), 
                "z": float(ligne[46:54].strip()) }

                dico_atoms[ligne[12:16].strip()] = dico_coord #Dictionnaire de dictionnaire 
                #Dictionnaire avec les atomes et pour chaque atome : dictionnaire de coordonnees
                dico_coord={}
                
                liste_atoms.append(dico_atoms)
                dico_atoms={}

    return liste_atoms

def calc_dist_3D(A,B):
    dst3d=math.sqrt((A[0]-B[0])**2+(A[1]-B[1])**2 +(A[2]-B[2])**2 )
    return dst3d



#La structure secondaire est définie par l'arrangement des liaisons hydrogène, 
#en conséquence, la définition exacte de celles-ci est cruciale. Dans DSSP, 
#la définition standard d'une liaison hydrogène dérive d'un modèle purement électrostatique. 
#DSSP attribue des charges partielles q1 de +0,42e et - 0,42e sur le carbone et l'oxygène du 
#carbonyle (C=O) et q2 de +0,20e et - 0,20e sur l'hydrogène et l'azote de l'amide (NH),


def calc_Energie(A,B):
    Q1_CO=0.42
    Q2_NH=-0.20
    




c=lecture_fichier_pdb("1cfc.pdb")
print(c)