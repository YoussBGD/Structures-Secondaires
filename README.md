


- Ce programme assigne des structures secondaires à des protéines en se basant sur les coordonnéeS des 
  atomes présentes dans des fichiers pdb avec la méthode DSSP.

- Ce programme fonctionne sur un seul fichier pdb (sur une seule frame ) ou sur une tajectoire de 
  dynamique moléculaire ( doit être sous format pdb).

- Notre projet est séparé en 4 fichiers python, ces fichiers se trouvent dans le dossier nommé "scripts":
  1- Un fichier nommé DSSP_functions contenant toutes les fonctions qui calculent et assignent les structures secondaire.

  2- Un fichier nommé Draw_ss.py qui contient une fonction permettant d'afficher sous forme de flèches les feuillets 
     beta et sous forme ressorts les hélices alpha d'une seule frame (que nous pouvons choisir).

  3- Un fichier nommé Draw_ss_trajectory qui contient une fonction permettant d'afficher les structures secondaires
     sous forme de bondes colorées selon le type de la structure secondaire pour toute une trajectoire. 

  4- Un fichier nommé DSSP_launcher qui est le fichier qui lance notre programme.

nous avons aussi 3 autres dossiers:
  - Un dossier nommé "PDB_file" qui doit contenir les fichiers .pdb à analyser.
    Le programme peut également télécharger directement les fichiers sur la pdb et les copier directement dans ce dossier.
  - Un dossier nommé "Output" qui contiendra tous les fichiers de sortie.
  - Un dossier nommé "conda_env"  qui permettra de créer l'environnement contenant les bibliothèques nécéssaires pour
    ce programme. 

Lancement du programme: 
1-  Création de l'environnement conda contenant les bibliothèques requises pour ce programme:
	- Ouvrir un terminal dans le dossier conda_env
	- Taper : conda env create -f DSSP_M2ISDD.yml
	- Puis :  conda activate DSSP_M2ISDD

2 - Ouvrir un terminal dans le dossier scripts.
3 - Taper la commande : "python DSSP_launcher -h" Pour voir toutes les options de ce programme.
   
    Les informations suivantes s'afficheront sur votre terminal:

            ##############################################################################################
	    #   Optional arguments:                                                                      #
	    # 	  -h, --help :   show this help message and exit                                         #
            #                                                                                            #
	    # 	  -m MODEL, --model MODEL : Model number (for RMN structure) - Default = 1               #  
	    #   	        	                                                                 # 
	    #	  -c CHAIN, --chain CHAIN : Chain number to analyze - Default = A                        #
	    #			                                                                         #
	    #	  -d, --downloadPDB :    Download file from PDB - Default = False                        #
            #                                                                                            #
	    #	  -nh, --no_addH     :       Use ambertools reduce to add H - Default = True             #
            #                                                                                            #
	    #	  -p, --plotSS   :       Secondary structure plot - Default = False                      # 
	    #			                                                                         #
	    #	  -trj, --DSSP_Trajectory : Do DSSP on a trajectory (on PDB file) - Default = False      #
	    #	                                                                                         #
            #     -stat,--statistiques :To calculate statistics of Secondary structures (on PDB file)    #
            #                           - Default = False		                                 #
	    #			                                                                         #
            ##############################################################################################

 
 Explication des différentes options, comment les utiliser et comment lancer notre programme:
---------------------------------------------------------------------------------------------

 A - (-d) 
     Premièrement il faut mettre le nom du fichier pdb que vous souhaitez étudier.
     Ici nous avons deux cas :
     - Premier cas, si le fichier pdb que vous souhaitez étudier se trouve déja dans le dossier 
     PDB_file, par exemple le fichier 1cfc.pdb, ici vous n'aurez qu'à taper la commande suivante :
       python DSSP_launcher 1cfc (+ la suite des options)
     - Deuxième cas, quand vous souhaitez étudier un fichier qui ne se trouve pas dans le dossier PDB_file, ici vous 
     pourez le télécharger facilement dans ce dossier en mettant l'option -d ou --downloadPDB après l'ID de votre 
     structure pdb que vous souhaitez étudier, la commande sera la suivante:
       python DSSP_launcher 1cfc -d (+ la suite des options)   
      
 
 B - (-m) 
      Pour choisir quel modèle vous souhaitez étudier dans un fichier pdb contenant plusieurs 
      modèles comme le fichier 1cfc.pdb ou un fichier trajectoire contenant plusieurs frames, vous 
      pouvez utiliser la commande -m ou --model. Si vous n'indiquez pas cette option le modèle 1 (ou la frame 1) 
      sera choisi par défaut.
      - Par exemple la ligne de commande pour choisir le modèle numéro 10 dans le fichier pdb 1cfc sera donc :
        python DSSP_launcher 1cfc -m 10 (+ la suite des options)

 C - (-c) 
     Pour choisir la chaine de la protéine dans le fichier pdb que vous souhaitez étudier, il faut utiliser 
     l'option -c ou --chain. Au cas ou vous ne préciseriez pas cette option c'est la chaine A qui sera choisi 
     par défaut (NB : Certains fichiers pdb ne posèdent pas de noms de chaines, ce programme est aussi adapté 
     pour ce cas et ce ne sera pas necessaire de préciser -c).
     - Par exemple si on souhaite étudier la chaine A du modèle 10 de notre fichier pdb 1cfc.pdb la commande sera 
       la suivante:
        python DSSP_launcher 1cfc -m 10 -c A (+ la suite des options)

 D - (-nH) 
     Notre programme ajoute les hydrogènes par défaut à la protéine du fichier pdb choisi avec Ambertools reduce. 
     Cependant si vous ne souhaitez pas ajouter d'hydrogènes vous pouvez mettre l'option -nh ou --no_addH (Il est conseillé de toujours laisser le
     programme protoner la protéine pour ne pas avoir d'anomalies lors de la recherche de structures secondaires). 

 E - (-trj) 
     Si vous souhaitez étudier les structures secondaires de toutes les frames d'un fichier pdb contenant une trajectoire de dynamique 
     moléculaire entière, il faut mettre l'option -trj ou --DSSP_Trajectory 
     Par exemple, pour etudier le fichier traj_bclxl.pdb contenant 50 frames, se trouvant dans le dossier PDB_file, il faut ecrire :
        python DSSP_launcher traj_bclxl -trj (+ la suite des options)

 F - (-p) 
     Pour avoir une image en sortie (soit pour une seule frame (un seul modèle) ou une trajectoire) afin de voir 
     les différentes structures secondaires de la protéine du fichier choisi, il faut mettre l'option -p ou --plotSS.
     - par exemple si vous souhaitez afficher les structures socondaire pour le fichier trajectoire "traj_bclxl" 
      il faut lancer la commande suivante : Python DSSP_launcher traj_bclxl -trj -p 
     - Si nous souhaitons le faire pour une seule frame (la 10 par exemple ) la commande sera comme la suivante:
        Python DSSP_launcher traj_bclxl -m 10 -p
 
 
 G - (-stat) 
     Si vous souhaitez avoir des statistiques sur le pourcentage de chaque structure secondaire dans la protéine étudiée,
     vous pouvez ajouter l'option -stat ou --statistiques. Un fichier nommé stat_ID-Pdb.txt sera crée
       python DSSP_launcher 1cfc -d -stat

 H - Par défaut, un fichier de sortie nommé SS_{ID-Pdb} contenant les types de structures secondaires que forment les résidus du fichier 
     étudié est crée (soit pour une trajectoire entière ou une seule frame (un seul modèle)). 
     Ce fichier contiendra le numéro du résidu et le type de structure secondaire qu'il forme, voici les types de
     structures secondaires que détecte notre programme:
     	
        E = (Extended) strand of a β-pleated sheet
	G = Right-handed 3_10 helix
	H = Right-handed α-helix
	I = Right-handed π-helix
	T = Turn
	C = Coil




	     
	       
     
