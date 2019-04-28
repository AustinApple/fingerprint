import pandas
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import openbabel
import numpy as np 
import os 
from biopandas.mol2 import PandasMol2
import openbabel

def smile2ECFP_2048_number(row):
    from rdkit.Chem import AllChem
    from rdkit import Chem
    from rdkit.Chem import Draw
    import pandas
    from collections import Counter
    m = Chem.MolFromSmiles(row['smile'])
    #print(row['smile'])
    info={}
    fp = AllChem.GetMorganFingerprintAsBitVect(m,2,nBits=2048,bitInfo=info).ToBitString()
    
    try:
        for bit in info:
            row["num_"+str(bit)] = len(info[bit])
        return row  
    except:
         print(row["smile"],"ECFP_2048 something is wrong!!")



def ECFP_2048(row):
    try:
        mol = Chem.MolFromSmiles(row['smile'])
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048).ToBitString()
        for i in range(2048):
            row[str(i)] = fp[i]
        return row  
    except:
         print(row["smile"],"ECFP_2048 something is wrong!!")

def SYBYL(row):
    try:
        obconversion = openbabel.OBConversion()
        obconversion.SetInAndOutFormats("smi", "mol2")
        mol = openbabel.OBMol()
        obconversion.ReadString(mol,row["smile"])  # read molecule from database 
        mol.AddHydrogens()
        output_mol2 = obconversion.WriteString(mol)  # transform smiles into mol2
        file = open("molecule.mol2","w+")   # write mol2 format into the file, molecule.mol2.
        file.write(output_mol2)
        file.close()
        molecule_mol2 = PandasMol2().read_mol2("molecule.mol2")  # use biopandas to static the discriptors
        for element in molecule_mol2.df['atom_type'].value_counts().index:
            if element == 'Al':
                row['Al'] = molecule_mol2.df['atom_type'].value_counts()['Al']
            if element == 'B':
                row['B'] = molecule_mol2.df['atom_type'].value_counts()['B']
            if element == 'Br':
                row['Br'] = molecule_mol2.df['atom_type'].value_counts()['Br']
            if element == 'C.1':
                row['C.1'] = molecule_mol2.df['atom_type'].value_counts()['C.1']
            if element == 'C.2':
                row['C.2'] = molecule_mol2.df['atom_type'].value_counts()['C.2']
            if element == 'C.3':
                row['C.3'] = molecule_mol2.df['atom_type'].value_counts()['C.3']
            if element == 'C.ar':
                row['C.ar'] = molecule_mol2.df['atom_type'].value_counts()['C.ar']
            if element == 'C.cat':
                row['C.cat'] = molecule_mol2.df['atom_type'].value_counts()['C.cat']
            if element == 'Ca':
                row['Ca'] = molecule_mol2.df['atom_type'].value_counts()['Ca']
            if element == 'Cl':
                row['Cl'] = molecule_mol2.df['atom_type'].value_counts()['Cl']
            if element == 'F':
                row['F'] = molecule_mol2.df['atom_type'].value_counts()['F']
            if element == 'H':
                row['H'] = molecule_mol2.df['atom_type'].value_counts()['H']  
            if element == 'Li':
                row['Li'] = molecule_mol2.df['atom_type'].value_counts()['Li']  
            if element == 'Mg':
                row['Mg'] = molecule_mol2.df['atom_type'].value_counts()['Mg']        
            if element == 'N.1':
                row['N.1'] = molecule_mol2.df['atom_type'].value_counts()['N.1']         
            if element == 'N.2':
                row['N.2'] = molecule_mol2.df['atom_type'].value_counts()['N.2']         
            if element == 'N.3':
                row['N.3'] = molecule_mol2.df['atom_type'].value_counts()['N.3']         
            if element == 'N.4':
                row['N.4'] = molecule_mol2.df['atom_type'].value_counts()['N.4']         
            if element == 'N.am':
                row['N.am'] = molecule_mol2.df['atom_type'].value_counts()['N.am']         
            if element == 'N.ar':
                row['N.ar'] = molecule_mol2.df['atom_type'].value_counts()['N.ar']        
            if element == 'N.pl3':
                row['N.pl3'] = molecule_mol2.df['atom_type'].value_counts()['N.pl3']         
            if element == 'Na':
                row['Na'] = molecule_mol2.df['atom_type'].value_counts()['Na']        
            if element == 'O.2':
                row['O.2'] = molecule_mol2.df['atom_type'].value_counts()['O.2']        
            if element == 'O.3':
                row['O.3'] = molecule_mol2.df['atom_type'].value_counts()['O.3'] 
            if element == 'O.co2':
                row['O.co2'] = molecule_mol2.df['atom_type'].value_counts()['O.co2']  
            if element == 'P.3':
                row['P.3'] = molecule_mol2.df['atom_type'].value_counts()['P.3']    
            if element == 'S.2':
                row['S.2'] = molecule_mol2.df['atom_type'].value_counts()['S.2']        
            if element == 'S.3':
                row['S.3'] = molecule_mol2.df['atom_type'].value_counts()['S.3']
            if element == 'S.O2':
                row['S.O2'] = molecule_mol2.df['atom_type'].value_counts()['S.O2']   
            if element == 'S.O':
                row['S.O'] = molecule_mol2.df['atom_type'].value_counts()['S.O']       
            if element == 'Si':
                row['Si'] = molecule_mol2.df['atom_type'].value_counts()['Si']    
            if element == 'Zn':
                row['Zn'] = molecule_mol2.df['atom_type'].value_counts()['Zn']       
        return row
    except:
        print(row["smile"],"SYBYL something is wrong!!")

def featurize(list_smile,do_smile2ECFP_2048_number=None, do_smile2_SYBYL_ECFP=None):
    if do_smile2ECFP_2048_number:
        smile_list_num_bit =['smile']+["num_"+str(i) for i in range(2048)]
        data = pd.DataFrame(columns=smile_list_num_bit)
        data['smile'] = list_smile
        data = data.apply(smile2ECFP_2048_number, axis=1).fillna(0)
        return data

    elif do_smile2_SYBYL_ECFP:
        smile_list_num_bit =['smile']+[str(i) for i in range(2048)]
        atom_type=['Al','B','Br','C.1','C.2','C.3','C.ar','C.cat','Ca','Cl','F','H',
                    'Mg','N.1','N.2','N.3','N.4','N.am','N.ar','N.pl3','Na',
                    'O.2','O.3','O.co2','P.3','S.2','S.3','S.O','S.O2','Si','Zn']
        columns_order = smile_list_num_bit+atom_type
        data = pandas.DataFrame(columns=columns_order)
        data['smile'] = list_smile
        data = data.apply(ECFP_2048, axis=1)        
        data = data.apply(SYBYL, axis=1).fillna(0)
        data = data[columns_order]
        return data


















    # print("I am here")
    # # arrange the columns 
    # folder = str(foldername)    
    # list_name= ["file_name","smile"]
    # list_num_bit =["num_"+str(i) for i in range(2048)]
    # # atom_type=['B','Br','C.1','C.2','C.3','C.ar','C.cat','Ca','Cl','F','H',
    # #              'Mg','N.1','N.2','N.3','N.4','N.am','N.ar','N.pl3','Na',
    # #              'O.2','O.3','P.3','S.2','S.3','S.O','S.O2','Si','Zn']
    # # columns_order = list_name+list_bit+atom_type
    # data = pandas.DataFrame(columns=list_name+list_num_bit)
    
    # # get the folder 
    # directory = os.getcwd()
    # absolute_folder_directory = os.path.join(directory,folder)
    # list_file = os.listdir(absolute_folder_directory)
    
    # list_smile = []
    # list_filename = []

    # for file in list_file:
    #     list_filename.append(file)
    #     obconversion = openbabel.OBConversion()
    #     obconversion.SetInAndOutFormats("mol", "smi")   #set input and output format 
    #     mol = openbabel.OBMol()
    #     obconversion.ReadFile(mol,os.path.join(absolute_folder_directory,file))  #input 
    #     smile = obconversion.WriteString(mol)[:-21]
    #     list_smile.append(smile)

    # data["smile"] = list_smile       # put the smiles into the dataframe 
    # data["file_name"] = list_filename  # put the file_name into the dataframe 
    # # featurization 
    # # data = data.apply(ECFP_2048, axis=1)        
    # # data = data.apply(SYBYL, axis=1).fillna(0)
    # data = data.apply(smile2ECFP_2048_number, axis=1).fillna(0)
    # data.to_csv(outputname)
