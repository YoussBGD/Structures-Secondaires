import math
import numpy as np
import pandas as pd


def lecture_fich_pdb(file_name):
    with open(file_name, "r") as f_pdb:
        coor_lst = []
        for ligne in f_pdb :
            if (ligne.startswith("ATOM") and ligne[12:16].strip() in ["H","H1","O","N","C"]) or (str(ligne[17:21].strip())=='PRO' and str(ligne[12:16].strip())=="HA"):
            
                # Création du dictionnaire vide.
                dico_infos = {}
    
                # Extraction de l'atome
                dico_infos["atom"] = str(ligne[12:16].strip())
                
                # Extraction du nom du résidu.
                dico_infos["residu"] = str(ligne[17:21])
                
                # Extraction du numéro du résidu.
                dico_infos["N° resid"] = int(ligne[22:26])
                
                # Extraction de la coordonnée x.
                dico_infos["X"] = float(ligne[30:38])
                
                # Extraction de la coordonnée y. 
                dico_infos["Y"] = float(ligne[38:46])
                
                # Extraction de la coordonnée z.
                dico_infos["Z"] = float(ligne[46:54])

                # Extraction de la chaine
                dico_infos["chaine"] = str(ligne[21:23])
                coor_lst.append(dico_infos)
                
        DF_coord = pd.DataFrame(coor_lst)
        
    return DF_coord
      


def dist_3D(A,B):
    dst3d=math.sqrt((A["X"]-B["X"])**2+(A["Y"]-B["Y"])**2 +(A["Z"]-B["Z"])**2 )
    return dst3d

def coord_un_res(data_frame,num_res):
    dic_coord={}
    dico_atoms={}
    liste_atom=[]
    if num_res==1:
        i=0
    else:
        i=4*(num_res-1) #pour trouver le numéro de ligne d'un atome on fait (numero du résidu-1)*4 
                  #(2 c'est les 3H du residu n ter -1 car on commence a 0 )
    while data_frame["N° resid"][i]==num_res:

        dic_coord={  #Dictionnaire de coordonnees
                "X": float(data_frame["X"][i]), 
                "Y": float(data_frame["Y"][i]), 
                "Z": float(data_frame["Z"][i])}

        dico_atoms[str(data_frame["atom"][i])] = dic_coord
        dic_coord={}
        i=i+1

    return dico_atoms

        
#n,c,o,h
#0 1 2 3 
#A: N et H
#B: C et O


#La structure secondaire est définie par l'arrangement des liaisons hydrogène, 
#en conséquence, la définition exacte de celles-ci est cruciale. Dans DSSP, 
#la définition standard d'une liaison hydrogène dérive d'un modèle purement électrostatique. 
#DSSP attribue des charges partielles q1 de +0,42e et - 0,42e sur le carbone et l'oxygène du 
#carbonyle (C=O) et q2 de +0,20e et - 0,20e sur l'hydrogène et l'azote de l'amide (NH),

def calc_Energie(A,B):
    Q1_CO=0.42
    Q2_NH=0.20
    if("H" in B.keys()):
        E=Q1_CO*Q2_NH*((1/dist_3D(B["N"],A["O"]))+(1/dist_3D(B["H"],A["C"]))-(1/dist_3D(B["H"],A["O"]))-(1/dist_3D(B["N"],A["C"])))*332
    else:
        if("H1" in B.keys()):
            E=Q1_CO*Q2_NH*((1/dist_3D(B["N"],A["O"]))+(1/dist_3D(B["H1"],A["C"]))-(1/dist_3D(B["H1"],A["O"]))-(1/dist_3D(B["N"],A["C"])))*332
        else:
            E=Q1_CO*Q2_NH*((1/dist_3D(B["N"],A["O"]))+(1/dist_3D(B["HA"],A["C"]))-(1/dist_3D(B["HA"],A["O"]))-(1/dist_3D(B["N"],A["C"])))*332
    return E

def calc_beta(data_frame,alpha):
    #truc=int(len(data_frame)/4)
    pont=[]
    list_beta=[]
    truc=141
    for x in range(1,truc):
        for y in range(1,truc):
            if([x,y] not in alpha and [y,x] not in alpha and abs(x-y)>10):
                A= coord_un_res(data_frame,x)
                B= coord_un_res(data_frame,y)
                E=calc_Energie(A,B)
                if(E<-0.5):
                    pont=[x,y,E]
                    list_beta.append(pont)
                pont=[]
    return list_beta


c=lecture_fich_pdb("1cfc.pdb")
#print(c["residu"][14780:14800])

print(c)
cc=coord_un_res(c,4)
print(cc)

print(len(c)/4)


alpha=[]
res=[]
for i in range(1,140):

    y=i+4
    atm1=coord_un_res(c,i)
    atm2=coord_un_res(c,y)
    E=calc_Energie(atm1,atm2)
    if E<-0.5:
        print('res1 : {} | res2 : {} | energie = {}'.format(i,y,E))
        res=[i,y]
        alpha.append(res)
        res=[]
print(alpha)


k=calc_beta(c,alpha)
print(k)
