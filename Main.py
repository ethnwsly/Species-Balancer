# Species Balancer
#
# Author: Ethan Chen
# California State University Fullerton
# Nikolaidis Lab
# Last Updated: May 31, 2018
# Python 3.5
# Required package(s): Biopython
#
# Purpose: This script takes in FASTA files and returns pruned FASTA files that contain only
# sequences that are from species common to EVERY FASTA file in your batch.
# Basically this insures that your FASTA files contain the same species across all files and removes species that
# that are not found in every file.
#
#
# How to use: Edit input folder file path and output folder file path below
#


import os
from Bio import SeqIO
from Gene import Gene


def findcommonset(mylist):

    if len(mylist) < 2:
        raise ValueError("Gene list is not large enough, < 2")

    answerset = mylist[0].getSpeciesList()

    for workinggene in mylist:
        workinglist = workinggene.getSpeciesList()
        answerset = [val for val in answerset if val in workinglist]
    return answerset


def main():

    # user_path = path for your Input folder
    # output_path = path for Output
    user_path = "/home/ethan/PycharmProjects/Species_Balancer/Input/"
    output_path = "/home/ethan/PycharmProjects/Species_Balancer/Output/"

    # Check if path given is valid, throws an exception if invalid
    assert os.path.exists(user_path)
    assert os.path.exists(output_path)

    # Change working directory to user path
    os.chdir(user_path)
    items = os.listdir(user_path)

    filelist = []
    for names in items:
        if names.endswith(".fa"):
            filelist.append(names)
        elif names.endswith(".fasta"):
            filelist.append(names)

    print("Files Imported: ")
    print(filelist)
    # filelist contains list of fasta file names

    genelist = []
    # genelist contains a list of Gene objects. Each Gene object contains a list of the species associated with the gene

    for gene in filelist: # gene = 1 fasta file
        records = list(SeqIO.parse(gene,"fasta")) # records = list of sequence entries
        #print("\n")
        #print("Found %i records" % len(records))
        species_list = []

        for indiv_record in records:
            # Parses species name from description
            # WARNING: Parsing will only work if the description format stays constant
            # If NCBI changes their description formatting, the parsing here will need to be changed as well

            species_name = indiv_record.description.split("[")[1][:-1]
            species_list.append(species_name)

        genelist.append(Gene(species_list))

    commonspeciesset = findcommonset(genelist)
    print("Common species found:")
    print(commonspeciesset)

    # Copy and edit each file in filelist to only contain records that fit the commonspeciesset
    # Output files in another folder

    for gene in filelist:
        records = list(SeqIO.parse(gene,"fasta"))
        new_records = []

        for indiv_record in records:
            species_name = indiv_record.description.split("[")[1][:-1]
            if species_name in commonspeciesset:
                new_records.append(indiv_record)

        os.chdir(output_path)
        SeqIO.write(new_records, gene, "fasta")
        os.chdir(user_path)

    print("Files written to output")


if __name__ == "__main__":
    main()















