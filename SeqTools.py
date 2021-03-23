#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 12:57:02 2020

@author: dal75
"""
import csv, os
import numpy as np
from datetime import datetime
from Bio import SeqIO


def read_covid_sequences(file_name):
    '''
    Reads the covid sequences metadata and makes a list of each entry as a dict
    '''
    seq_names = []
    with open(file_name, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            seq_names.append(row)
        for i in range(len(seq_names)):
            seq_names[i]['sample_date'] = datetime.strptime(seq_names[i]['sample_date'],"%Y-%m-%d")
            seq_names[i]['epi_week'] = int(seq_names[i]['epi_week'])
            try:
                seq_names[i]['lineage_support'] = float(seq_names[i]['lineage_support'])
            except:
                pass
            
    return seq_names

def sequence_selector(seq_names,**kwargs):
    """
    Return new list of sequences that match filter
    Filters:
        country = UK-ENG, UK-WLS, UK-SCT, UK-NIR
        startdate = datetime object
        enddate = datetime object
        area = relevant only to England and Scotland
        epi_week = list of epiweeks to consider (either int, list of int, or np.array of int)
        lineage = lineage code
    """
    country = kwargs.get("country",False)
    startdate = kwargs.get("startdate",False)
    enddate = kwargs.get("enddate",False)
    area = kwargs.get("area",False)
    epi_week = kwargs.get("epi_week",False)
    lineage = kwargs.get("lineage",False)
    
    if (startdate and enddate):
        #check that if there are both a startdate and enddate, that they are correct
        #chronilogically and also both datetimes
        try:
            if (startdate > enddate):
                print("Error: startdate must be before enddate")
                return
        except:
            print("Either startdate or enddate are not of datetime format")
            return
            
    if(country):
        #check if country os correct format and filters for country
        ctrs = ["UK-ENG","UK-WLS","UK-NIR","UK-SCT"]
        if (country in ctrs):
            new_files = []
            for f in seq_names:
                if f['adm1'] == country:
                    new_files.append(f)
            seq_names = new_files
        else:
            print("Incorrect country format, formats are {0}, {1}, {2}, {3}".format(ctrs[0],ctrs[1],ctrs[2],ctrs[3]))
    
    if(startdate):
        #check for datetime format and then filter for all seqs after startdate
        if type(startdate) == datetime:
            new_files = []
            for f in seq_names:
                if f['sample_date'] >= startdate:
                    new_files.append(f)
            seq_names = new_files
        else:
            print("startdate is not of datetime format")
            
    if(enddate):
        #check for datetime format and then filter for all seqs before enddate
        if type(enddate) == datetime:
            new_files = []
            for f in seq_names:
                if f['sample_date'] <= enddate:
                    new_files.append(f)
            seq_names = new_files
        else:
            print("startdate is not of datetime format")
    
    if (epi_week):
        if type(epi_week) == int:
            weeks = np.array([epi_week], dtype = int)
        elif type(epi_week) == list:
            try:
                weeks = np.array(epi_week, dtype = int)
            except:
                print("List does not containe solely ints")
                return
        elif type(epi_week) == np.ndarray:
            if (epi_week.dtype != int):
                print("np.array is not all ints")
                return
            else:
                weeks = epi_week
        else:
            print("epi_week is of incorrect type")
            return
        new_files = []
        for f in seq_names:
            if f['epi_week'] in weeks:
                new_files.append(f)
        seq_names = new_files
                
            
        
        
    return seq_names

def get_areas(files):
    '''
    From files gives list of areas
    '''
    
    areas = {"England":[],"Scotland":[],"Wales":[],"Northern_Ireland":[]}
    
    for f in files:
        name = f['sequence_name'].split("/")
        if (name[0] == "Scotland"):
            ar = name[1][:3]
            if not any(char.isdigit() for char in ar):
                if ar not in areas[name[0]]:
                    areas[name[0]].append(name[1][:3])
        else:
            ar = name[1][:4]
            if not any(char.isdigit() for char in ar):
                if ar not in areas[name[0]]:
                    areas[name[0]].append(name[1][:4])
                
    new_areas = {}
    new_areas['UK-ENG'] = areas['England']
    new_areas['UK-SCT'] = areas['Scotland']
    new_areas['UK-WLS'] = areas['Wales']
    new_areas['UK-NIR'] = areas['Northern_Ireland']
    return new_areas
        
def fetch_fastas(filtered_seqs, cog_file, folder_name):
    '''
    
    Parameters
    ----------
    filtered_seqs : list of dicts
        details on sequence names that have been filteres
    cog_file : fasta file location
        location of fasta file to extract data from
    folder_name : string
        what to name the folder of the output fastas

    Returns
    -------
    None.

    '''
    print("starting process")
    #get list of sequence names to select
    names = []
    for seq in filtered_seqs:
        names.append(seq['sequence_name'])
    
    #add matching sequences to a list
    new_seqs = []
    with open(cog_file, newline='') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.name in names:
                new_seqs.append(record)
    print("Completed reading in sequences")
    #checks/creates folder to store fastas
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name, 0o777)
    if not os.path.isdir(os.path.join(folder_name,"genomes")):
        os.mkdir(os.path.join(folder_name,"genomes"), 0o777)
            
        print("Creating folders for storage")
        
    #creats a metadata file for the fastas
    with open(os.path.join(folder_name,"metadata.csv"),'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=filtered_seqs[0].keys())
        writer.writeheader()
        for seq in filtered_seqs:
            writer.writerow(seq)
    with open(os.path.join(folder_name,"sequences.txt"),'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='	')
        writer.writerow(["seq_name","aln_name","seq_path","annotation_path"])
        i = 0
        for seq in new_seqs:
            seq_name = "_".join(seq.name.split("/"))
            aln_name = "seq"+str(i)
            seq_path = os.path.join(os.getcwd(),folder_name,"genomes",seq_name+".fasta")
            annotation_path = "NA"
            writer.writerow([seq_name,aln_name,seq_path,annotation_path])
            i += 1
            print("Writing sequence "+str(i)+" to sequences.txt")
    for seq in new_seqs:
        with open(os.path.join(folder_name,"genomes","_".join(seq.name.split("/"))+".fasta"), 'w',newline='') as file:
            SeqIO.write(seq, file, 'fasta')
    print("Finished writing separate fasta files")
        
    return

def load_fastas(folder_name):
    return
    


