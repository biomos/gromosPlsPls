"""
This file is part of GROMOS.

Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
See <https://www.gromos.net> for details.

GROMOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

##############################################################################
# buildHybrid: A Python tool for Generating Hybrid Building Blocks with the 
# GROMOS force field.
#
# A documentation can be found in the doc-project: doc/misc/Docu_buildHybrid/
#
# Authors: Bartosz Stankiewicz, Niels Hansen (hansen@itt.uni-stuttgart.de)
#          University of Stuttgart, Germany
#
# buildHybrid is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details:
# <http://www.gnu.org/licenses/>.
##############################################################################
#
# Modifications : [list below changes after May 2017 - with initials, location 
# and date - and update time stamp]
#
# NH, 05/31/2017: initial version
#
###############################################################################
#!/usr/bin/python
from subprocess import call
from subprocess import check_call
import sys, re, time, copy, fileinput, os.path

# Generates a hybrid .mtb file out of two residues and perturbation topologies
def main(show):

    # Specify the respective locations:  
    forceFileName = (
     "...your-path-to.../forcefields/official/54a7.mtb")
    massFileName = (
     "...your-path-to.../forcefields/official/54a7.ifp")

    run(show,forceFileName,massFileName)

def showRes():
    # Aromatic residues: 15,16,17,18,19,20,26,27,30,31.  
    # At the end of each residue line you'll find the one/two letter code.  
    print ("\nThis version is intended for beta testing. Please check carefully\n"
           + "all files that were generated and report any errors to the "
           + "authors.")
    print "\n---Building a hybrid residue---\n"
    print "List of available residues to choose from:"
    print "0=Alanine\t\t\t\t\t\tALA\tA"
    print "1=Arginine (protonated; charge +e)\t\t\tARG\tR"
    print "2=Arginine (deprotonated; neutral)\t\t\tARGN\tRN"
    print "3=Asparagine\t\t\t\t\t\tASN\tN"
    print "4=Asparagine (coordinated with ZN)\t\t\tASN1\tN1"
    print "5=Aspartic acid (deprotonated; charge -e)\t\tASP\tD"
    print "6=Aspartic acid (protonated; neutral)\t\t\tASPH\tDH"
    print "7=Cysteine (deprotonated; charge -0.5e)\t\t\tCYS\tC"
    print "8=Cysteine (protonated; neutral)\t\t\tCYSH\tCH"
    print "9=Cysteine (1st member of S-S bridge)\t\t\tCYS1\tC1"
    print "10=Cysteine (2nd member of S-S bridge)\t\t\tCYS2\tC2"
    print "11=Glutamine\t\t\t\t\t\tGLN\tQ"
    print "12=Glutamine acid (deprotonated; charge -e)\t\tGLU\tE"
    print "13=Glutamine acid (protonated; neutral)\t\t\tGLUH\tEH"
    print "14=Glycine\t\t\t\t\t\tGLY\tG"
    print "15=Histidine (protonated at ND1; neutral)\t\tHISA\tHA"
    print "16=Histidine (protonated at NE2; neutral)\t\tHISB\tHB"
    print "17=Histidine (protonated at ND1 and NE2; charge +e)\tHISH\tHH"
    print "18=Histidine (coupled to HEME at NE2; neutral)\t\tHIS1\tH1"
    print "19=Histidine (coupled to HEMC at NE2: neutral)\t\tHIS2\tH2"
    print "20=Hydroxyproline\t\t\t\t\tHYP\tPH"
    print "21=Isoleucine\t\t\t\t\t\tILE\tI\n22=Leucine\t\t\t\t\t\tLEU\tL"
    print "23=Lysine (deprotonated; neutral)\t\t\tLYS\tK"
    print "24=Lysine (protonated; charge +e)\t\t\tLYSH\tKH"
    print "25=Methionine\t\t\t\t\t\tMET\tM"
    print "26=Phenylalanine\t\t\t\t\tPHE\tF\n27=Proline\t\t\t\t\t\tPRO\tP"
    print "28=Serine\t\t\t\t\t\tSER\tS\n29=Threonine\t\t\t\t\t\tTHR\tT"
    print "30=Tryptophan\t\t\t\t\t\tTRP\tW"
    print "31=Tyrosine\t\t\t\t\t\tTYR\tY\n32=Valine\t\t\t\t\t\tVAL\tV\n"  

def getRes():
    print("Please define a hybrid residue, i.e. type in the number for "
              + "first and second residue:")
    while True:
        try:
                 resOne = int(raw_input("1st residue: "))
                 break
        except ValueError:
            print "Not a number! Try again."
    while True:
        try:
                 resTwo = int(raw_input("2nd residue: "))
                 break
        except ValueError:	
            print "Not a number! Try again."
    return [resOne,resTwo]

def run(show,forceFileName,massFileName):
    [resOne, resTwo, verbose] = processArguments(show)
    [withGly, withAla, withPro, hypPro, withArg] = setCases(resOne, resTwo)

    renamedAtoms = []
    rememberAI = []
    exclusSho = []
    atomInfoHyb = []
    exclusHyb = []
    exclusHybUp = []
    excludedAtoms = []
    newBondInt = []
    oldAmount = -1
    newBonds = []
    if not withPro:
        # Initializes all needed files and lists.
        argFileName = "add_atom.arg"
        [forceFile,argFile,residues,aromatic] = initialize(
                                              forceFileName, argFileName,
                                              withArg)
        # Determines the the side chain length and chooses longer one as
        # basis.  
        [longerRes,lastNumSho,lastNumLo,forceLines] = basisCheck(
            resOne, resTwo, forceFile, residues, aromatic)
        # Writing an add_atom.arg file.  
        if withArg:
                writeArg(
                 longerRes, lastNumSho, argFile,
                 forceFileName, residues, withArg, withPro, 4, 4)
                # .arg has to be closed in order to execute add_atom.  
                argFile.close()
        # Writing to .mtb by executing add_atom from this python program and
        # writing and manipulating output file.
        if longerRes==resOne:
            shorterRes=resTwo
        else:
            shorterRes=resOne 
        outFileName = output(longerRes, shorterRes, resOne, resTwo,
                             forceLines, residues, argFileName,
                             withArg, withPro, False)
        # Determining atom information and exclusion list of shorter
        # residue.  
        print "Residues: "+residues[resOne][1]+" and "+residues[resTwo][1]
        print "Chosen residue basis: "+residues[longerRes][1]
        [atomInfo,exclusSho,bondsSho] = excluSho(
            shorterRes, forceLines, residues, lastNumSho, withArg, hypPro)

        # Comparing atom names and returning changed atom names, if equal.  
        [replace,contentsMtb] = compareAtomN(
            outFileName, atomInfo, residues, shorterRes, withPro)
        if not withArg:
            contentsMtb.insert(-1,"#  IB   JB  MCB  NCO  IND  CON\n")
        # Determining exclusion list of hybrid
        [atomInfoHyb,exclusHyb,bondsLo,bondsStr] = excluLo(
            longerRes, contentsMtb, residues, withPro, hypPro)
        rememberAI = copy.deepcopy(atomInfo)
        rememberBo = copy.deepcopy(bondsLo)
        # Number of lengths, etc, before updating
        oldAmount = copy.deepcopy(bondsLo[0])
        # Changing the atom names in shorter residue
        renamedAtoms = changeAtomN(replace,atomInfo)

        # Updating the exclusion list of the hybrid with the exclusions in
        # shorter residue.  
        if not hypPro:
            exclusHybUp = updateExclu(exclusSho, exclusHyb, 2,
                                      True, [], hypPro)

            # Adding atoms of hybrid to exclusion list of atoms in short 
            # residue.  
            excludedAtoms = updateAtom(atomInfo, lastNumSho,
                                       lastNumLo, hypPro)
        # Returning an integer list of lists with all added lengths, angles, 
        # etc.  
        newBondInt = updateBonds(atomInfo, bondsSho, bondsLo,
                                 withArg, withPro)
        # Making a string from the updated bonds
        newBonds = bondsToString(
            newBondInt, bondsStr, residues, shorterRes, withArg)

        # Changing RNME, putting together all lines from initial.mtb and all
        # generated strings and saving modified.mtb.  
        atomsToFile(
            exclusHybUp, excludedAtoms, newBonds, contentsMtb, residues,
            shorterRes, longerRes, outFileName, lastNumSho, replace, withArg,
            hypPro)

        # Writing a perturbation file for each direction
        case = perturb(
            massFileName, residues, forceLines, shorterRes, longerRes,
            lastNumSho, lastNumLo, exclusSho, exclusHyb, rememberAI,
            atomInfoHyb, renamedAtoms, withArg, withPro, hypPro, withAla,
            withGly)

        # Printing for testing
        printAll(
        renamedAtoms, False, rememberAI, False, exclusSho, False,
        atomInfoHyb, False, exclusHyb, False, exclusHybUp, False,
        excludedAtoms, False, newBondInt, False, 3, oldAmount,
        newBonds, False)
    else:
        # Initializes all needed files and lists.
        argFileName1 = "add_atom_1.arg"
        argFileName2 = "add_atom_2.arg"
        [forceFile, argFile1, argFile2, residues, aromatic] = initializePro(
            forceFileName, argFileName1, argFileName2)
        [longerRes, lastNumSho, lastNumLo, forceLines] = basisCheckPro(
            resOne, resTwo, forceFile, residues, aromatic, withArg)
        # add_atom_1.arg
        writeArg(longerRes, lastNumSho, argFile1,
                 forceFileName, residues, withArg, withPro, 1, 1)
        argFile1.close()
        if longerRes==resOne:
            shorterRes=resTwo
        else:
            shorterRes=resOne
        outFileName1 = output(longerRes, shorterRes, resOne, resTwo,
                              forceLines, residues, argFileName1,
                              withArg, withPro, False)
        outFileManip = manipArgOne(
                        outFileName1, longerRes, shorterRes, residues)
        # add_atom_2.arg.  Last argument of writeArg corresponds to the atom
        # after which 'New' atoms will be added.  This is 6 in PRO and 8 in
        # HYP, because there was already inserted H.  Second last argument
        # corresponds to number of atoms in shorter residue, that will NOT
        # be part of hybrid.  
        if not withGly and longerRes == 20:
            # Basis is HYP
            if withArg:
                writeArg(longerRes, lastNumSho, argFile2,
                         outFileManip, residues, withArg, withPro, 4, 8)
            else:
                writeArg(longerRes, lastNumSho, argFile2,
                         outFileManip, residues, withArg, withPro, 3, 8)
        elif not withGly and longerRes == 27:
            # Basis is PRO
            if withArg:
                writeArg(longerRes, lastNumSho, argFile2,
                         outFileManip, residues, withArg, withPro, 4, 6)
            else:
                writeArg(longerRes, lastNumSho, argFile2,
                         outFileManip, residues, withArg, withPro, 3, 6)
        else:
            # GLY <-> HYP, PRO
            argFileName2 = argFileName1
        argFile2.close()
        outFileName2 = output(longerRes, shorterRes, resOne, resTwo,
                              forceLines, residues, argFileName2,
                              withArg, withPro, True)
        # Determining atom information and exclusion list of shorter
        # residue.  
        print "Residues: "+residues[resOne][1]+" and "+residues[resTwo][1]
        print "Chosen residue basis: "+residues[longerRes][1]
        [exclusSho, exclusHyb, atomInfoHyb, renamedsAtoms] = manipInitPro(
           longerRes, shorterRes, lastNumSho, 
           lastNumLo, forceLines, residues, 
           withArg, withPro, outFileName2,
           hypPro, withGly)

        # Writing a perturbation file for each direction
        case = perturb(
            massFileName, residues, forceLines, shorterRes, longerRes,
            lastNumSho, lastNumLo, exclusSho, exclusHyb, rememberAI,
            atomInfoHyb, renamedAtoms, withArg, withPro, hypPro,
            withAla, withGly)
    forceFile.close()
    # If user doesn't want to keep processing files
    if not verbose:
        delete(shorterRes, longerRes, residues)
    # Updates the final .mtb file with a header
    header(shorterRes, longerRes, residues, case)
    # Truncates the final .mtb to X_Y.mtb
    rename(shorterRes, longerRes, residues, case)
    # Uncomment for printing the generated files:
    # 
    # printFiles(shorterRes, longerRes, residues, 
    #           True, 70, True,
    #           False, case, False)

# Handles arguments and returns verbose, resOne and resTwo.  
def processArguments(show):
    verbose = False
    numbers = False
    lenShow = len(show)
    red = 0
    try:
        check = int(show[-2])
        check = int(show[-1])
        numbers = True
        red = -2
    except (ValueError, IndexError) as e:
        try:
            check = int(show[-1])
            sys.exit("Invalid arguments or number of arguments. Exiting.")
        except (ValueError, IndexError) as e:
            # There are no residues in argument list
            numbers = False
            red = -2
    iCount = 1
    # Check, if all non-residue arguments are valid
    while iCount < lenShow+red:
        if show[iCount] != "-s" and show[iCount] != "-v": 
            sys.exit("Invalid arguments or number of arguments. Exiting.")
        iCount = iCount + 1
    
    if not "-s" in show:
        showRes()
    if "-v" in show:
        verbose = True
    if not numbers:
        [resOne, resTwo] = getRes()
    else:
        resOne = int(show[-1])
        resTwo = int(show[-2])
    return [resOne, resTwo, verbose]

def setCases(resOne, resTwo):
    # Program assumes no Glycine, Alanine or (Hydroxy)Proline by default.   
    withGly = False
    withAla = False
    withPro = False
    hypPro = False
    if (resOne == 27 and resTwo == 20) or (resOne == 20 and resTwo == 27):
        hypPro = True
    else:
        if resOne == 27 or resTwo == 27 or resOne == 20 or resTwo == 20:
            withPro = True
            if resOne==27 or resTwo== 27:
                print ("\nWARNING: A perturbation of bonds in the pyrrolidine "
                   + "ring is not yet\nimplemented. Therefore the non-proline "
                   + "endstate may likely\nhave artificial conformational "
                   + "restrictions.\n")
            else:
                print ("\nWARNING: A perturbation of bonds in the pyrrolidine "
                   + "ring is not yet\nimplemented. Therefore the "
                   + "non-hydroxyproline endstate may likely\nhave artificial "
                   + "conformational restrictions.\n")
        if resOne == 0 or resTwo == 0:
            withAla = True
        if resOne == 14 or resTwo == 14:
            withGly = True

    # Program writes .arg-File(s) by default
    withArg = True
    if hypPro == True:
        withArg = False
    else:
        if withGly and not (withPro or withAla):
            withArg = False
        elif withAla and not (withPro or withGly):
            withArg = False
        elif withAla and withGly and not withPro:
            withArg = False
    return [withGly, withAla, withPro, hypPro, withArg] 

# Initializes an array with both the one(two) and more letter code of all the
# residues.  Opens file for reading and writing.  
def initialize(forceFileName,argFileName,withArg):
    forceFile = open(forceFileName, 'r')
    if withArg:
        argFile = open(argFileName, 'w')
    else:
        argFile    = ""
    [residues,aromatic] = getResAro()
    return [forceFile,argFile,residues,aromatic]

# (Hydroxy) Proline requires another initializing procedure, because there are
# two argument files.  It also differs by the parameter withArg, because there
# needs to be executed add_atom regardless of whether Glycine is part of the
# hybrid or not.   
def initializePro(forceFileName, argFileName1, argFileName2):
    forceFile = open(forceFileName, 'r')
    argFile1 = open(argFileName1, 'w')
    argFile2 = open(argFileName2, 'w')
    [residues,aromatic] = getResAro()
    return [forceFile, argFile1, argFile2, residues, aromatic]


def getResAro():
    residues = [
                 ['A', 'ALA', 'Alanine'],
                 ['R','ARG','Arginine (protonated; charge +e)'],
                 ['RN','ARGN','Arginine (deprotonated; neutral)'],
                 ['N','ASN','Asparagine'],
                 ['N1','ASN1','Asparagine (coordinated with ZN)'],
                 ['D','ASP','Aspartic acid (deprotonated; charge -e)'],
                 ['DH','ASPH','Aspartic acid (protonated; neutral)'],
                 ['C','CYS','Cysteine (deprotonated; charge -0.5e)'],
                 ['CH','CYSH','Cysteine (protonated; neutral)'],
                 ['C1','CYS1','Cysteine (1st member of S-S bridge)'],
                 ['C2','CYS2','Cysteine (2nd member of S-S bridge)'],
                 ['Q','GLN','Glutamine'],
                 ['E','GLU','Glutamine acid (deprotonated; charge -e)'],
                 ['EH','GLUH','Glutamine acid (protonated; neutral)'],
                 ['G','GLY','Glycine'],
                 ['HA','HISA','Histidine (protonated at ND1; neutral)'],
                 ['HB','HISB','Histidine (protonated at NE2; neutral)'],
                 ['HH','HISH',
                       'Histidine (protonated at ND1 and NE2; charge +e)'],
                 ['H1','HIS1','Histidine (coupled to HEME at NE2; neutral)'],
                 ['H2','HIS2','Histidine (coupled to HEMC at NE2: neutral)'],
                 ['PH','HYP','Hydroxyproline'], ['I','ILE','Isoleucine'],
                 ['L','LEU','Leucine'],
                 ['K','LYS','Lysine (deprotonated; neutral)'],
                 ['KH','LYSH','Lysine (protonated; charge +e)'],
                 ['M','MET','Methionine'], ['F','PHE','Phenylalanine'],
                 ['P','PRO','Proline'], ['S','SER','Serine'], 
                 ['T','THR','Threonine'],
                 ['W','TRP','Tryptophan'], ['Y','TYR','Tyrosine'],
                 ['V','VAL','Valine']
                 ]
    aromatic = [15,16,17,18,19,20,26,27,30,31]
    return [residues,aromatic]

# Compares the length of two residues and returns natural number
# representation of the longer residue as well as number of the last atom
# inside the shorter residue.  
# 
# Does also a check of the residue numbers.  Returns also the read out force
# file.  
def basisCheck(resOne, resTwo, forceFile, residues, aromatic):
    # Checks, if given numbers are valid
    if resOne<0 or resTwo<0 or resOne>32 or resTwo>32 or resOne==resTwo:
        sys.exit("Defined residues numbers invalid or equal. Exiting.")
    # Searches for information block of the residues in forceFile and
    # determines last atom number before trailing atoms.  
    forceLines = forceFile.readlines()
    numLines=len(forceLines)
    chainNum1=-1
    chainNum2=-1
    for iCount in range(numLines):
        if residues[resOne][1] == forceLines[iCount].rstrip():
            trailing="# trailing atoms".rstrip()
            for jCount in range(numLines):
                if trailing == forceLines[iCount+jCount].rstrip():
                    words=forceLines[iCount+jCount-1].split()
                    chainNum1=int(words[0])
                    break
        if chainNum1>-1:
            break
    # Does the same for second residue
    for iCount in range(numLines):
        if residues[resTwo][1] == forceLines[iCount].rstrip():
            trailing="# trailing atoms".rstrip()
            for jCount in range(numLines):
                if trailing == forceLines[iCount+jCount].rstrip():
                    words=forceLines[iCount+jCount-1].split()
                    chainNum2=int(words[0])
                    break
        if chainNum2>-1:
            break
    # Checking and comparing side chain lengths
    if chainNum1==-1:
        sys.exit("First chain side not determined. Exiting.")
    if chainNum2==-1:
        sys.exit("Second chain side not determined. Exiting.")
    if chainNum1>chainNum2:
        return [resOne,chainNum2,chainNum1,forceLines]
    elif chainNum1<chainNum2:
        return [resTwo,chainNum1,chainNum2,forceLines]
    else:
        # Equal side chain lengths. take the aromatic one
        if resOne in aromatic:
            return [resOne,chainNum2,chainNum1,forceLines]
        elif resTwo in aromatic:
            return [resTwo,chainNum1,chainNum2,forceLines]
        else:
            return [resOne,chainNum2,chainNum1,forceLines]

# In case of (Hydroxy) Proline, returns:   
# [NNR long, last #atom short, last #atom long, force lines].  
# NNR stands for 'natural number representation' of the residue.  
def basisCheckPro(resOne, resTwo, forceFile, residues, aromatic, withArg):
    if resOne<0 or resTwo<0 or resOne>32 or resTwo>32 or resOne==resTwo:
        sys.exit("Defined residues numbers invalid or equal. Exiting.")
    # Searches for information block of the residues in forceFile and
    # determines last atom number before trailing atoms.  
    forceLines = forceFile.readlines()
    numLines = len(forceLines)
    numChainLong = -1
    numChainShort = -1
    searchRes = -1
    nnrl = -1
    if resOne == 20 and not resTwo == 27:
        # Basis is Hydroxyproline.  resTwo is regular residue.  
        nnrl = 20
        numChainLong = 7
        searchRes = resTwo
    elif resTwo == 20 and not resOne == 27:
        # Basis is Proline.  resOne is regular residue.  
        nnrl = 20
        numChainLong = 7
        searchRes = resOne
    elif resOne == 27 and not resTwo == 20:
        # Basis is Proline.  resTwo is regular residue.  
        nnrl = 27
        numChainLong = 5
        searchRes = resTwo
    elif resTwo == 27 and not resOne == 20:
        # Basis is Hydroxyproline.  resOne is regular residue.  
        nnrl = 27
        numChainLong = 5
        searchRes = resOne
    else:
        # Either residue 1 or residue 2 is Proline, whereas the other one must
        # be Hydroxyproline.
        nnrl = 20
        numChainLong = 7
        numChainShort = 5
    if searchRes >= 0:
        if withArg:
            for iCount in range(numLines):
                if residues[searchRes][1] == forceLines[iCount].rstrip():
                    trailing="# trailing atoms".rstrip()
                    for jCount in range(numLines):
                       if trailing == forceLines[iCount+jCount].rstrip():
                            words = forceLines[iCount+jCount-1].split()
                            numChainShort = int(words[0])
                            break
                if numChainShort>-1:
                    break
        else:
            # Short residue is Glycine
            numChainShort = 3
    return [nnrl, numChainShort, numChainLong, forceLines]

# Writes an .arg input-file with all parameters for add_atom and the hybrid.  
# 
# Careful: lastNum is not the last atom number of residues[numRes], but of the
# shorter residue.  backbone is the number of atoms up to the carbon beta
# atom in numRes. If PRO or HYP, it is the last side chain atom number.  
def writeArg(
        numRes, lastNum, argFile, forceFileName,
        residues, withArg, withPro, backbone, start):
    if withPro and withArg:
        pro = 1
    else:
        pro = 0
    argFile.write("@file "+forceFileName+"\n")
    resName    = residues[numRes][1]
    argFile.write("@build "+resName+"\n")
    argFile.write("@start "+str(start)+"\n")
    addNumAtom = lastNum-backbone+pro
    if backbone > 1:
            argFile.write("@number "+str(addNumAtom)+"\n\n")
    else:
            # Backbone <= 1 indicates .arg file for first add_atom execution
            argFile.write("@number "+str(1)+"\n\n")

# Writes the initial .mtb file and returns shorter residue natural number
# and .mtb file name.  
def output(
        longerRes, shorterRes, resOne, resTwo, forceLines,
        residues, argFileName, withArg, withPro, second):
    if not withPro:
        outFileName = (residues[shorterRes][1] + "_" + residues[longerRes][1]
                           + "_initial.mtb")
    else:
        if not second:
            outFileName = (residues[shorterRes][1]
                           + "_" + residues[longerRes][1]
                           + "_initial_1.mtb")
        else:
            outFileName = (residues[shorterRes][1]
                           + "_" + residues[longerRes][1]
                           + "_initial_2.mtb")
    if withArg or withPro:
        outFile    = open(outFileName, 'w')    
        check_call("add_atom @f "+ argFileName, shell=True, stdout=outFile)
        outFile.close()
    else:
        outFile    = open(outFileName, 'w')    
        contentsMtb = getBuildingBlock(forceLines,residues,longerRes)
        for iCount in contentsMtb:
            outFile.write(iCount),
        outFile.close()
    return outFileName

def getBuildingBlock(forceLines, residues, res):
    # iCount points to residue identifier, e.g. PHE
    iCount = 0
    # jCount points to END of building block
    jCount = 0
    while residues[res][1] != forceLines[iCount].rstrip():
        iCount = iCount + 1
    jCount = iCount
    while "END" != forceLines[jCount].rstrip():
        jCount = jCount + 1
    contentsMtb = []
    contentsMtb.append(forceLines[iCount-4])
    contentsMtb.append(forceLines[iCount-2])
    contentsMtb.append(forceLines[iCount-1])
    contentsMtb.append(forceLines[iCount])
    while iCount != jCount:
        iCount = iCount + 1
        if "#@FREELINE\n" != forceLines[iCount]:
            contentsMtb.append(forceLines[iCount])
    return contentsMtb

# Adds a TITLE block to a .mtb file
def manipArgOne(outFileName, longerRes, shorterRes, residues):
    mtbFile = open(outFileName, 'r')
    contentsMtb = mtbFile.readlines()
    mtbFile.close()
    outFileManip = (residues[shorterRes][1] + "_" + residues[longerRes][1]
                           + "_modified_1.mtb")
    mtbFile = open(outFileManip, 'w')
    contentsMtb.insert(0,"TITLE\n")
    iCount = 1
    while contentsMtb[iCount].split()[0]=="#":
        iCount = iCount + 1
    contentsMtb.insert(iCount, "END\n")
    for jCount in contentsMtb:
        mtbFile.write(jCount)
    mtbFile.close()
    return outFileManip

# This functions is only called, if HYP or PRO are involved
def manipInitPro(longerRes, shorterRes, lastNumSho, lastNumLo,
                 forceLines, residues, withArg, withPro, 
                 outFileName, hypPro, withGly):
    outFile = open(outFileName, 'r')
    contents = outFile.readlines()
    outFile.close()
    # 1) Head of the hybrid file
    head = contents[0:13]
    direction = residues[shorterRes][0] + "_" + residues[longerRes][0] + "\n"
    # 2) Upper body of hybrid file
    body1 = contents[14:19]
    prec = []
    # This condition is equivalent of asking: is HYP part of hybrid?  
    if longerRes == 20:
        prec.append("   -1                                 "
                        + "5    0    1    2    3    8\n")
    # This condition is only true, if PRO is the only residue in hybrid.  
    else:
        prec.append("   -1                                 "
                        + "5    0    1    2    3    6\n")
    # 3) prec holds the preceding exclusion list lines
    prec.append(contents[20])
    # 4) Middle body of the file 
    body2 = contents[21:23]

    [atomInfo, exclusSho, bondsSho] = excluSho(shorterRes, forceLines,
                                                 residues, lastNumSho,
                                                 withArg, hypPro)
    # Leading atom numbers must be corrected in case of HYP or PRO
    jCount = 0
    if longerRes == 20:
        add = 4
    else:
        add = 2
    oldNumbers = []
    while jCount < len(atomInfo):
        if len(atomInfo[jCount].split())>6:
            temp = atomInfo[jCount][6:]
            oldNumber = atomInfo[jCount].split()[0]
            oldNumbers.append(int(oldNumber))
            newNumber = (getEmpty(5,len(str(int(oldNumber)+add)))
                         + str(int(oldNumber)+add) + " ")
            atomInfo[jCount]= newNumber + temp
        jCount = jCount + 1
    # Comparing atom names and returning changed atom names, if equal.  
    [replace,contentsMtb] = compareAtomN(
            outFileName, atomInfo, residues, shorterRes, withPro)
    # Determining exclusion list of hybrid.  exclusHyb consists in case of
    # PRO of N, H, CA, CB, CG and CD (6 atoms) and in case of HYP of N, H,
    # CA, CB, CG, OD1, HD1 and CD2 (8 atoms)
    [atomInfoHyb,exclusHyb,bondsLo,bondsStr] = excluLo(
            longerRes, contentsMtb, residues, withPro, hypPro)
    # Ignores the Hydrogen atom
    del exclusHyb[1]
    numbers = getExclInt(exclusHyb,2)

    # Hardcoded N-parameters with already updated exclusion list.  
    # 5) enn = N lines, eidsh = H lines.  Careful: The Nitrogen gets charge of
    # shorter residue.  
    # newNum = 1 N / 2 H / ... / 4 CB / # 5 ..... 10 # = 6 
    newNum = lastNumSho - 4
    if longerRes == 20:
        firstNew = 9
        addedExcluH = range(3,firstNew)
        addedExcluH.insert(0,len(addedExcluH))
    else:
        firstNew = 7
        addedExcluH = range(3,firstNew)
        addedExcluH.insert(0,len(addedExcluH)) 
    # PRO: rangeNew = 7, 8, ... , 15.  Add these exclusions to N.   
    rangeNew = range(firstNew,firstNew + newNum + 2)
    # Add H-atom number
    rangeNew.insert(0,2)
    addedExcluN = addExclNum(rangeNew,numbers[0])
    temp = addedExcluN[1:]
    temp.sort()
    temp.insert(0,addedExcluN[0])
    addedExcluN = copy.deepcopy(temp)
    blockSho = getBuildingBlock(forceLines, residues, shorterRes)
    lenBlock = len(blockSho)
    found = False
    jCount = 0
    while jCount < lenBlock and not found:
        if "    1 N" == blockSho[jCount][0:7]:
            found = True
        else:
            jCount = jCount + 1
    splitN = blockSho[jCount].split()
    found = False
    jCount = 0
    while jCount < lenBlock and not found:
        if "    2 H" == blockSho[jCount][0:7]:
            found = True
        else:
            jCount = jCount + 1
    splitH = blockSho[jCount].split()
    enn = getAtom(
                  int(splitN[0]), splitN[1], int(splitN[2]),
                  int(splitN[3]), splitN[4], int(splitN[5]),
                  addedExcluN)
    eidsh = getAtom(
                    int(splitH[0]), splitH[1], int(splitH[2]),
                    int(splitH[3]), splitH[4], int(splitH[5]),
                    addedExcluH)


    # This is the atom number list for all backbone atom exclusions to add.  
    # newAtomsNum is the list 7, 8, 9, ... of new atoms.  
    newAtomsNum = range(firstNew, firstNew + newNum + 1)
    if withArg:
        rangeNew = range(firstNew,firstNew + newNum + 1)
    else:
        rangeNew = range(firstNew,firstNew + newNum)
    addedExclu = rangeNew

    # 6) ca holds CA lines and cb holds CB lines
    caLine = []
    found = False
    iCount = 0
    while iCount < len(exclusHyb) and not found:
        if "    3 CA" == exclusHyb[iCount][0:8]:
            found = True
            caLine.append(exclusHyb[iCount])
        else:
            iCount = iCount + 1
    while iCount != len(exclusHyb)-1 and len(exclusHyb[iCount+1].split()) < 7:
        caLine.append(exclusHyb[iCount+1])
    splitCA = exclusHyb[iCount].split()

    cbLine = []
    found = False
    iCount = 0
    while iCount < len(exclusHyb) and not found:
        if "    4 CB" == exclusHyb[iCount][0:8]:
            found = True
            cbLine.append(exclusHyb[iCount])
        else:
            iCount = iCount + 1
    while iCount != len(exclusHyb)-1 and len(exclusHyb[iCount+1].split()) < 7:
        cbLine.append(exclusHyb[iCount+1])
        iCount = iCount + 1
    splitCB = exclusHyb[iCount].split()
    numbersCA = getExclInt(caLine,1)
    numbersCA = numbersCA[0]
    numbersCB = getExclInt(cbLine,1)
    numbersCB = numbersCB[0]
    ca = getAtom(
                    int(splitCA[0]), splitCA[1], int(splitCA[2]),
                    int(splitCA[3]), splitCA[4], int(splitCA[5]),
                    addExclNum(addedExclu,numbersCA))
    cb = getAtom(
                    int(splitCB[0]), splitCB[1], int(splitCB[2]),
                    int(splitCB[3]), splitCB[4], int(splitCB[5]),
                    addExclNum(addedExclu,numbersCB))
    exclusHybUp = updateExclu(exclusSho, exclusHyb[1:3], 2, True, [], hypPro)

    # 7) ceeGee holds CG line(s)
    cgLine = []
    iCount = 0
    found = False
    while iCount < len(exclusHyb) and not found:
        if "    5 CG" == exclusHyb[iCount][0:8]:
            found = True
            cgLine.append(exclusHyb[iCount])
        else:
            iCount = iCount + 1
    while iCount != len(exclusHyb)-1 and len(exclusHyb[iCount+1].split()) < 7:
        cgLine.append(exclusHyb[iCount+1])
    numbersCG = getExclInt(cgLine,1)
    numbersCG = numbersCG[0]
    splitH = exclusHyb[iCount].split()
    ceeGee = getAtom(
                    int(splitH[0]), splitH[1], int(splitH[2]),
                    int(splitH[3]), splitH[4], int(splitH[5]),
                    addExclNum(addedExclu,numbersCG))
    cd = []
    if longerRes == 20:
        hypProAtom = ["    6 OD1","    7 HD1","    8 CD2"]
    else:
        hypProAtom = ["    6 CD "]
    hypProLen = len(hypProAtom)
    for jCount in range(hypProLen):
        temp = []
        iCount = 0
        found = False
        while iCount < len(exclusHyb) and not found:
            if hypProAtom[jCount] == exclusHyb[iCount][0:9]:
                found = True
                temp.append(exclusHyb[iCount])
            else:
                iCount = iCount + 1
        while (iCount < len(exclusHyb)-1
                and len(exclusHyb[iCount+1].split()) < 7):
            temp.append(exclusHyb[iCount+1])
        numbersTemp = getExclInt(temp,1)
        numbersTemp = numbersTemp[0]
        splitTemp = exclusHyb[iCount].split()
        tempAtom = getAtom(
            int(splitTemp[0]), splitTemp[1], int(splitTemp[2]),
            int(splitTemp[3]), splitTemp[4], int(splitTemp[5]),
            addExclNum(addedExclu,numbersTemp))
        cd.append(tempAtom)

    iCount = 0
    atomInfo = []
    if withArg and not withGly:
        while"    4 CB" != blockSho[iCount][0:8]:
            iCount = iCount + 1
        while "# trailing atoms" != blockSho[iCount].rstrip()[0:16]:
            atomInfo.append(blockSho[iCount])
            iCount = iCount + 1
    rememberAI = copy.deepcopy(atomInfo)
    atomList = []
    atomList.append(enn)
    atomList.append(eidsh)
    atomList.append(ca)
    atomList.append(cb)
    atomList.append(ceeGee)
    for iCount in range(len(cd)):
        atomList.append(cd[iCount])
    numbersNew = getExclInt(atomInfo, len(newAtomsNum))
    if longerRes == 20:
        shift = 5
    else:
        shift = 3
    for iCount in range(len(numbersNew)):
        for jCount in range(1,len(numbersNew[iCount])):
            numbersNew[iCount][jCount] = numbersNew[iCount][jCount] + shift
    newAtomList = []
    iCount = 0
    jCount = 0
    while iCount < len(newAtomsNum):
        while len(atomInfo[jCount].split()) < 7:
            jCount = jCount + 1
        splitTemp = atomInfo[jCount].split()
        tempAtom = getAtom(
                int(splitTemp[0])+shift, splitTemp[1], int(splitTemp[2]),
                int(splitTemp[3]), splitTemp[4], int(splitTemp[5]),
                numbersNew[iCount])
        newAtomList.append(tempAtom)
        iCount = iCount + 1
        jCount = jCount + 1
    [replace,contentsMtb] = compareAtomN(
            outFileName, newAtomList, residues, shorterRes, withPro)
    renamedAtoms = changeAtomN(replace,newAtomList)

    # Shifting bondsSho.  Do not touch bondsSho[0] and bondsSho[...][-1].  
    # The problem is, that every bond/angle/dihedral, that contains Nitrogen
    # (2) or the sidechain starting from atom CA (3), must be added.  bondsSho
    # must be updated, though, i.e. shifted, because the function compares 
    # bondsSho to newAtomList, where the new atoms are appended to the end of
    # the side chain of the hybrid.
    for iCount in bondsSho:
        iCount,
    for iCount in range(1,len(bondsSho)):
        for kCount in range(len(bondsSho[iCount])):
            tooSmall = True
            lenBond = len(bondsSho[iCount][kCount])-1
            jCount = 0
            while tooSmall and jCount < lenBond:
                if bondsSho[iCount][kCount][jCount] < 2:
                    jCount = jCount + 1
                else:
                    tooSmall = False
            if not tooSmall:
                for jCount in range(len(bondsSho[iCount][kCount])-1):
                    if bondsSho[iCount][kCount][jCount] > 3:
                        bondsSho[iCount][kCount][jCount] = (
                         bondsSho[iCount][kCount][jCount] + shift)
    excludedAtoms = []
    oldAmount = copy.deepcopy(bondsLo[0])
    newBondInt = updateBonds(newAtomList, bondsSho, bondsLo, withArg, withPro)
    newBonds = bondsToString(newBondInt, bondsStr, residues, shorterRes,
                                 withArg)
    jCount = 0
    mtbLen = len(contentsMtb)
    while jCount < mtbLen and contentsMtb[jCount] != "# trailing atoms\n":
        jCount = jCount + 1
    iCount = 0
    while iCount < mtbLen and contentsMtb[iCount] != "# bonds\n":
        iCount = iCount + 1
    body3 = contentsMtb[jCount:iCount]
    sizeRep = len(replace)
    sizeAtomList = len(atomList)
    iCount = 0
    for iCount in range(sizeRep):
        # Searches for shifted atom entry
        atomNum = (getEmpty(5,len(str(4 + int(replace[iCount][0]))))
                + str(4 + int(replace[iCount][0])) + " ")
        jCount = 0
        found = False
        while jCount < sizeAtomList and not found:
            if len(atomList[jCount].split()) < 7:
                jCount = jCount + 1
            else:
                if atomList[jCount].split()[1] == replace[iCount][1][0:-1]:
                    found = True
                else:
                    jCount = jCount + 1
        replAtom = replace[iCount][1][0:-1]+"A"
        sub = replAtom + getEmpty(7, len(replAtom))
        atomList[jCount] = atomList[jCount][0:6] + sub + atomList[jCount][13:]
        iCount = iCount + 1
    atomsToFile = []
    for iCount in head:
        atomsToFile.append(iCount)
    for iCount in direction:
        atomsToFile.append(iCount)
    for iCount in body1:
        atomsToFile.append(iCount)
    for iCount in prec:
        atomsToFile.append(iCount)
    for iCount in body2:
        atomsToFile.append(iCount)
    for iCount in atomList:
        atomsToFile.append(iCount)
    atomsToFile.append("### " + residues[shorterRes][1] + " atoms\n")
    for iCount in newAtomList:
        atomsToFile.append(iCount)
    atomsToFile.append("###\n")
    for iCount in body3:
        atomsToFile.append(iCount)
    for iCount in newBonds:
        atomsToFile.append(iCount)
    atomsToFile.append("# LJ exceptions\n# NEX\n    0\n"+
     "#  IB   JB  MCB  NCO  IND  CON\nEND")
    outFile = open((residues[shorterRes][1] + "_" + residues[longerRes][1]
            + "_modified_2.mtb"),'w')
    for iCount in range(len(atomsToFile)):
        outFile.write(atomsToFile[iCount])
    # Printing for testing
    printAll(
        renamedAtoms, False, rememberAI, False, exclusSho, False,
        atomInfoHyb, False, exclusHyb, False, exclusHybUp, False,
        excludedAtoms, False, newBondInt, False, 1, oldAmount,
        newBonds, False)
    return [exclusSho, exclusHyb, atomInfoHyb, renamedAtoms]

# Gets atoms of a building block from fifth atom and returns atomInfo.  
# Gets exclusion list of carbon alpha and beta of this building block and
# returns exclus.  
def getAtomInfo(forceLines, residues, res, withArg, hypPro):
    numLinesForce = len(forceLines)
    # Searching for the atoms of the res building block and inserting the
    # whole line (breaking lines inclusively) into one list element in list
    # atomInfo.  Copies the exclusion entries of CA and CB of shorter residue
    # in list exclus.  
    atomInfo=[]
    exclus=[]
    leave=False
    exLineN=1
    for iCount in range(numLinesForce):
        if residues[res][1] == forceLines[iCount].rstrip():
            # Searching for the 4 CB entry and 3 CA entry. withArg is True,
            # if hybrid does not consist of Glycine.  
            if withArg or res == 0:
                if not hypPro:
                    trailingCB="CB"
                else:
                    trailingCB="CG"
                kCount=0
                jCount=0
                lCount=0
                for jCount in range(numLinesForce):
                    splitted=forceLines[iCount+jCount].split()
                    if len(splitted)>1:
                        carbonBeta=splitted[1]
                    else:
                        carbonBeta=None
                    if trailingCB == carbonBeta:
                        # Found 4 CB entry. remember the line in the force
                        # field, where the atom information of the shorter
                        # residue begins.  
                        jCount=iCount+jCount+1
                        break
                for kCount in range(numLinesForce):
                    # Leave the loop, if trailing atom entry is reached, else
                    # copy all lines.  
                    if ("# trailing atoms".rstrip()
                            == forceLines[jCount+kCount].rstrip()):
                        leave=True
                        # Copy the last line before trailing atom entry and
                        # then leave.  
                        atomInfo.append(forceLines[jCount+kCount-1])
                        break
                    atomInfo.append(forceLines[jCount+kCount-1])
                # Searches upwards for carbon alpha
                for kCount in range(numLinesForce):
                    trailingCA="CA"
                    splitted=forceLines[jCount-kCount-1].split()
                    if len(splitted)>1:
                        carbonAlpha=splitted[1]
                    else:
                        carbonAlpha=None
                    exclus.insert(0,forceLines[jCount-kCount-1])
                    # Checking, if CB entry has an exclusion list consisting
                    # of more than one line and extracting these.  
                    # 
                    # Assumption: if the lines consist of less than 7 entries,
                    # it must be a broken line of exclusions, because a "real"
                    # atom line consists of at least 7 entries.  Furthermore,
                    # ensuring that CB is not the last of the side chain 
                    # atoms.  exLineN saves how many additional lines the
                    # exclusion list has.  
                    if kCount==0:
                        for mCount in range(numLinesForce):
                            splitted=forceLines[jCount+mCount].split()
                            if (len(splitted)<7 and 
                                "# trailing atoms".rstrip()
                                != forceLines[jCount+mCount].rstrip()):
                                exclus.append(forceLines[jCount+mCount])
                                exLineN=exLineN+1
                            else:
                                break
                    if trailingCA == carbonAlpha:
                        # Found 3 CA entry. Remember the line in the force
                        # field, where the atom information of the shorter 
                        # residue begins.  
                        lCount=jCount-kCount
                        break
            # Glycine is one of the hybrid residues.  atomInfo is
            # empty and exclusion list consists of only one atom,
            # the carbon alpha atom.  Adding a second line, that
            # pretends to be a CB line.  withArg might be also False,
            # if ALA is involved.  In this case, do not add a dummy
            # line.  
            else:
                trailingCA="CA"
                kCount=0
                jCount=0
                lCount=0
                for jCount in range(numLinesForce):
                    splitted=forceLines[iCount+jCount].split()
                    if len(splitted)>1:
                        carbonAlpha=splitted[1]
                    else:
                        carbonAlpha=None
                    if trailingCA == carbonAlpha:
                        # Found 4 CA entry. Remembering the line of 'CA atom'
                        # line  
                        jCount=iCount+jCount
                        break
                for kCount in range(numLinesForce):
                    # Leave the loop, if trailing atom entry is reached, else
                    # copy all lines.  
                    if ("# trailing atoms".rstrip()
                        == forceLines[jCount+kCount].rstrip()):
                        leave=True
                        # Copy the last line before trailing atom entry and
                        # then leave
                        break
                    exclus.insert(0,forceLines[jCount+kCount])
                exclus.insert(1,"    4 CB     15    4    0.00000   0   0\n")
                atomInfo.append("")
        if leave is True:
            break
    # Delete CB (or CD, if hypPro) entry in atomInfo
    delAtom=[]
    if not hypPro:
        search = 'CB'
    else:
        search = 'CG'
    for iCount in range(exLineN):
        splitted=atomInfo[iCount].split()
        if len(splitted) > 1 and splitted[0]=='4' and splitted[1]==search:
            for jCount in range(exLineN):
                del atomInfo[0]
    return [atomInfo,exclus]

# Returns carbon alpha, carbon beta, the side chain and the additional bonds
# of shorterRes.  
def excluSho(shorterRes, forceLines, residues, lastNumSho, withArg, hypPro):
    missing=range(5,lastNumSho+1)
    numLinesForce = len(forceLines)
    # Read in the force field file
    [atomInfo,exclus] = getAtomInfo(forceLines, residues, shorterRes,
                                    withArg, hypPro)
    # Removing all entries, that are not newly added to hybrid and decreasing
    # amount correspondingly.  
    for iCount in range(2):
        splitIt=exclus[iCount].split()
        exclus[iCount]=exclus[iCount][0:36]
        length=len(splitIt)
        lost=0
        for jCount in range(7,length):
            atom = int(splitIt[jCount-lost])
            if atom not in missing:
                del splitIt[jCount-lost]
                lost=lost+1
        splitIt[6]=str(int(splitIt[6])-lost)
        for jCount in range(6,len(splitIt)):
            empty=" "
            spaces=len(splitIt[jCount])
            if jCount==6:
                for kCount in range(2-spaces):
                    empty=empty+" "
            else:
                for kCount in range(4-spaces):
                    empty=empty+" "
            exclus[iCount]=exclus[iCount]+empty+splitIt[jCount]
        exclus[iCount]=exclus[iCount]+"\n"
    # Searching for bonds and angles
    bonds=[[], [], [], [], []]
    searchList=[]
    searchList.append("# bonds")
    searchList.append("# bond angles")
    searchList.append("# improper dihedrals")
    searchList.append("# dihedrals")
    searchList.append("#@FREELINE")
    # There are 4 blocks, and there will be a first column, where the amount
    # of each block is saved (index 0).  
    for mCount in range(1,5):
        leave=False
        for iCount in range(numLinesForce):
            if residues[shorterRes][1] == forceLines[iCount].rstrip():
                # Searching for the 'number of bonds' line
                searchIn=searchList[mCount-1]
                kCount=0
                jCount=0
                lCount=0
                for jCount in range(numLinesForce):
                    splitted=forceLines[iCount+jCount].split()
                    if len(splitted)>2:
                        searchTerm=splitted[0]+" "+splitted[1]+" "+splitted[2]
                    elif len(splitted)>1:
                        searchTerm=splitted[0]+" "+splitted[1]
                    else:
                        searchTerm=None
                    if searchIn == searchTerm:
                        # Found 'number of bonds' line. Remembering the line
                        # in the force field, where it reads the number of
                        # bonds.  
                        jCount=iCount+jCount+2
                        bonds[0].append(int(forceLines[jCount]))
                        break
                for kCount in range(3,numLinesForce):
                    # Leave the loop, if 'number of degrees of freedom' line
                    # entry is reached, else copy all lines except the second
                    # (##  IB   JB  MCB...).  
                    if (searchList[mCount].rstrip()
                        == forceLines[jCount+kCount].rstrip()):
                        leave=True
                        # Copy the last line before bond angles entry and then
                        # leave.  
                        bonds[mCount].append(forceLines[jCount+kCount-1])
                        break
                    bonds[mCount].append(forceLines[jCount+kCount-1])
            if leave is True:
                break
    # Transforming the bond strings into integers
    for lCount in range(1,len(bonds)):
        for kCount in range(len(bonds[lCount])):
            bonds[lCount][kCount]=[
                int(mCount) for mCount in bonds[lCount][kCount].split()
                ]
                
    # Inserting comments to distinguish shorter residue atoms
    comment="### "+residues[shorterRes][1]+" atoms\n"
    return [atomInfo,exclus,bonds]

# Determines exclusion entries of CA and CB of longer residue in list exclus
# from initial.mtb.  Determines side chain of longer residue.  
def excluLo(longerRes,contentsMtb,residues, withPro, hypPro):
    numLinesMtb = len(contentsMtb)
    atomInfo=[]
    exclus=[]
    leave=False
    exLineN=1
    trailingCB = ""
    begin = ""
    for iCount in range(numLinesMtb):
        if residues[longerRes][1] == contentsMtb[iCount].rstrip():
            if withPro:
                if longerRes == 20:
                    trailingCB = "CD2"
                    begin = '8'
                else:
                    trailingCB = "CD"
                    begin = '6'
            elif hypPro:
                trailingCB = "CG"
                begin = '5'
            else:
                # Searching for the '4 CB' entry and 3 CA entry
                trailingCB = "CB"
                begin = '4'
            kCount=0
            jCount=0
            lCount=0
            for jCount in range(numLinesMtb):
                splitted=contentsMtb[iCount+jCount].split()
                if len(splitted)>1:
                    carbonBeta=splitted[1]
                else:
                    carbonBeta=None
                if trailingCB == carbonBeta:
                    # Found 4 CB entry. remember the line in the hybrid file,
                    # where the atom information of the shorter residue begins
                    jCount=iCount+jCount+1
                    break
            for kCount in range(numLinesMtb):
                # Leave the loop, if trailing atom entry is reached, else copy
                # all lines
                if ("# trailing atoms".rstrip()
                    == contentsMtb[jCount+kCount].rstrip()):
                    leave=True
                    # Copy the last line before trailing atom entry and then
                    # leave.  
                    atomInfo.append(contentsMtb[jCount+kCount-1])
                    break
                atomInfo.append(contentsMtb[jCount+kCount-1])
            for kCount in range(numLinesMtb):
                if withPro or hypPro:
                    trailingCA = "N"
                else:
                    trailingCA = "CA"
                splitted=contentsMtb[jCount-kCount-1].split()
                if len(splitted)>1:
                    carbonAlpha=splitted[1]
                else:
                    carbonAlpha=None
                exclus.insert(0,contentsMtb[jCount-kCount-1])
                # Checking, if CB entry has an exclusion list consisting of
                # more than one line and extracting these. 
                # 
                # Assumption: if the lines consist of less than 7 entries,
                # it must be a broken line of exclusions,
                # because a "real" atom line consists of at least 7 entries.  
                # Furthermore, ensuring that CB is not the last of the side
                # chain atoms.  exLineN saves how many additional lines the
                # exclusion list has.  
                if kCount==0:
                    for mCount in range(numLinesMtb):
                        splitted=contentsMtb[jCount+mCount].split()
                        if (len(splitted)<7 and "# trailing atoms".rstrip()
                                != contentsMtb[jCount+mCount].rstrip()):
                            exclus.append(contentsMtb[jCount+mCount])
                            exLineN=exLineN+1
                        else:
                            break
                if trailingCA == carbonAlpha:
                    # Found 3 CA entry. Remembering the line in the force
                    # field, where the atom information of the shorter residue
                    # begins.  
                    lCount=jCount-kCount
                    break
        if leave is True:
            break
    # Searching for bonds and angles
    bonds=[[], [], [], [], []]
    searchList=[]
    searchList.append("# bonds")
    searchList.append("# bond angles")
    searchList.append("# improper dihedrals")
    searchList.append("# dihedrals")
    searchList.append("# LJ exceptions")
    # There are 4 blocks, and there will be a first column, where the amount
    # of each block is saved (index 0).
    for mCount in range(1,5):
        leave=False
        for iCount in range(numLinesMtb):
            if residues[longerRes][1] == contentsMtb[iCount].rstrip():
                # Searching for the 'number of bonds' line
                searchIn=searchList[mCount-1]
                kCount=0
                jCount=0
                lCount=0
                for jCount in range(numLinesMtb):
                    splitted=contentsMtb[iCount+jCount].split()
                    if len(splitted)>2:
                        searchTerm=splitted[0]+" "+splitted[1]+" "+splitted[2]
                    elif len(splitted)>1:
                        searchTerm=splitted[0]+" "+splitted[1]
                    else:
                        searchTerm=None
                    if searchIn == searchTerm:
                        # Found 'number of bonds' line. Remembering the line
                        # in the force field, where it reads the 'number of
                        # bonds' line.  
                        jCount=iCount+jCount+2
                        bonds[0].append(int(contentsMtb[jCount]))
                        break
                for kCount in range(3,numLinesMtb):
                    # Leave the loop, if 'number of degrees of freedom' line
                    # entry is reached, else copy all lines except the second
                    # (##  IB   JB  MCB...).  
                    if (searchList[mCount].rstrip()
                        == contentsMtb[jCount+kCount].rstrip()):
                        leave=True
                        # Copy the last line before bond angles entry and then
                        # leave.  
                        bonds[mCount].append(contentsMtb[jCount+kCount-1])
                        break
                    bonds[mCount].append(contentsMtb[jCount+kCount-1])
            if leave is True:
                break
    # Make string copy of bonds before transforming
    bondsStr = copy.deepcopy(bonds[1:])
    # Transforming the bond strings into integers
    for lCount in range(1,len(bonds)):
        for kCount in range(len(bonds[lCount])):
            bonds[lCount][kCount]=[
                         int(mCount) for mCount 
                         in bonds[lCount][kCount].split()
                         ]
    # Delete CB entry in atomInfo
    delAtom=[]
    for iCount in range(exLineN):
        splitted=atomInfo[iCount].split()
        if splitted[0]==begin and splitted[1]==trailingCB:
            for jCount in range(exLineN):
                del atomInfo[0]
    return [atomInfo,exclus,bonds,bondsStr]

# Determines atom list of longer residue and compares to insert atom names
# with those atom names.  Returns a list of all names to be changed and the 
# relative position.  Reads initial.mtb.  
def compareAtomN(outFileName, atomInfo, residues, shorterRes, withPro):
    # Make a complete copy of atomInfo to stay save
    copyAtom = atomInfo[:]
    try:
        outFile    = open(outFileName, 'r+')
        contentsMtb = outFile.readlines()
        outFile.close()
    except IOError:
        sys.exit("No initial.mtb file.")
    # Read in the .mtb file to be changed and closing the file
    numCon=len(contentsMtb)
    jCount=0
    atoms=[]
    # Find first atom and save line in jCount
    for iCount in range(numCon):
        if ("#ATOM ANM  IACM MASS        CGMICGM MAE MSAE".rstrip()
                == contentsMtb[iCount].rstrip()):
            jCount = iCount+1
            break
    for iCount in range(numCon):
        if "# trailing atoms".rstrip() == contentsMtb[jCount+iCount].rstrip():
            break
        else:
        # Ignoring all lines, that don't contain alphabet characters,
        # i.e. all exclusion list broken lines.  
            match=re.search('[a-zA-Z]', contentsMtb[jCount+iCount].rstrip())
            if match:
                splitted = contentsMtb[jCount+iCount].split()
                atoms.append(splitted[1])
    # Reducing copyAtom to just the atom names of hybrid input
    for iCount in range(len(copyAtom)):
        copyAtom[iCount] = copyAtom[iCount].split()
        if len(copyAtom[iCount])>1:
            copyAtom[iCount] = copyAtom[iCount][1]

    # Doing the actual comparison and renaming
    replace=[]
    for iCount in range(len(atoms)):
        for jCount in range (len(copyAtom)):
            if atoms[iCount]==copyAtom[jCount]:
                # If jCount==1 then copyAtom[jCount] is atom 6 xxx   yy.  
                replace.append([jCount,copyAtom[jCount]+"B"])
    return [replace,contentsMtb]

# Replaces the clashing atom names with additional residue letter codes.  
def changeAtomN(replace,atomInfo):
    replaceLen=len(replace)
    for iCount in range(len(replace)):
        spaces=" "
        spaceLen = 7-len(replace[iCount][1])
        for jCount in range(spaceLen-1):
            spaces=spaces+" "
        # Ignoring comment line (replace[iCount][0]+1)
        atomInfo[replace[iCount][0]] = (atomInfo[replace[iCount][0]][0:6]
                + replace[iCount][1] + spaces
                + atomInfo[replace[iCount][0]][13:])
    print "Note: "+str(replaceLen)+" atom name(s) changed"
    return atomInfo

# Given a list 'exclus' of 'atomAm' atoms, returns an integer list of all
# exclusions with the number of exclusions at index 0.  
def getExclInt(exclus, atomAm):
    numbers=[]
    posLine=0
    for jCount in range(atomAm):
        entries=[]
        entered=False
        for iCount in range(posLine,len(exclus)):
            lineSho=exclus[iCount].split()
            if len(lineSho)>6:
                if entered:
                    posLine=iCount
                    break
                entered=True
                # Copying from column 6
                for kCount in range(6,len(lineSho)):
                    entries.append(int(lineSho[kCount]))
            else:
                for kCount in range(len(lineSho)):
                    entries.append(int(lineSho[kCount]))
        numbers.append(entries)
    return numbers

# Given a preceding number, an atom name, an IAC, a mass, a charge, a charge
# member integer and an integer exclusion list, returns an atom string
# (multi)line.  
def getAtom(num, name, iac, mass, charge, cm, numbers):
    exclNum = numbers[0]
    atomString = (getEmpty(5,len(str(num))) + str(num) + " " 
                  + name + getEmpty(5,len(name)) + " " 
                  + getEmpty(3,len(str(iac))) + str(iac) + " "
                  + getEmpty(4,len(str(mass))) + str(mass) + " "
                  + getEmpty(10,len(charge)) + charge + " "
                  + getEmpty(3,len(str(cm))) + str(cm) + " "
                  + getEmpty(3,len(str(exclNum))) + str(exclNum))
    iCount = 1
    while iCount <= exclNum:
        atomString = (atomString + getEmpty(5,len(str(numbers[iCount])))
              + str(numbers[iCount]))
        if iCount%6 == 0 and iCount != exclNum:
            atomString = (atomString 
                + "\n                                       ")
        iCount = iCount + 1
    atomString = atomString + "\n"
    return atomString

# Unifies addList with excl, updates excl[0] by new elements.  
def addExclNum(add,base):
    addCopy = copy.deepcopy(add)
    baseCopy = copy.deepcopy(base[1:])
    addCopy = eraseDuplicates(addCopy)
    baseCopy = eraseDuplicates(baseCopy)
    result = addCopy + baseCopy
    result.sort()
    result = eraseDuplicates(result)
    result.insert(0,len(result))
    return result

def eraseDuplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

# Returns the updated exclusion list of the hybrid. atomAm stores the amount
# of atoms.  If no second exclusion list is given, it will be extracted.  If
# noNum=False, then exclusHyb must be empty list.  
def updateExclu(exclusSho, exclusHyb, atomAm, noNum, numbersHyb, hypPro):
    numbersSho = getExclInt(exclusSho,atomAm)
    posLine=0
    if noNum:
        if hypPro:
            numbersHyb = getExclInt(exclusHyb,4)
        else:
            numbersHyb = getExclInt(exclusHyb,atomAm)
    # Comparing the number lists
    insertPos=0
    # Running over all atoms to be added
    for kCount in range(atomAm):
        numbersHyb[kCount] = addExclNum(
                numbersSho[kCount][1:],numbersHyb[kCount])
    if noNum and hypPro:
        additional = [5]
        for iCount in range(len(numbersHyb)):
            numbersHyb.insert(iCount, addExclNum(additional, numbersHyb[iCount]))
            del numbersHyb[iCount+1]
    # Building the new exclusion list
    exclusHybUp = []
    if noNum:
        add=0
        for iCount in range(len(numbersHyb)):
            while len(exclusHyb[iCount+add].split())<7:
                add=add+1
            if numbersHyb[iCount][0]>6:
                exclusHybUp.append(exclusHyb[iCount+add][0:35])
                for jCount in range(len(numbersHyb[iCount])):
                    empty=""
                    lenDigit=len(str(numbersHyb[iCount][jCount]))
                    # Determining how many spaces must be put in list
                    if jCount==0:
                        spaces = 4-lenDigit
                    else:
                        spaces = 5-lenDigit
                    for kCount in range(spaces):
                        empty=empty+" "
                    exclusHybUp[iCount] = (exclusHybUp[iCount] + empty
                         + str(numbersHyb[iCount][jCount]))
                exclusHybUp[iCount]=exclusHybUp[iCount]+"\n"
            else:
                exclusHybUp.append(exclusHyb[iCount+add][0:35])
                for jCount in range(len(numbersHyb[iCount])):
                    empty=""
                    lenDigit=len(str(numbersHyb[iCount][jCount]))
                    # Determining how many spaces must be put in list
                    spaces = 5-lenDigit
                    if jCount==0:
                        spaces = 4-lenDigit
                    else:
                        spaces = 5-lenDigit
                    for kCount in range(spaces):
                        empty = empty+" "
                    exclusHybUp[iCount] = (exclusHybUp[iCount] + empty
                         + str(numbersHyb[iCount][jCount]))
                exclusHybUp[iCount]=exclusHybUp[iCount]+"\n"
    else:
        add=0
        for iCount in range(len(numbersHyb)):
            while len(exclusSho[iCount+add].split())<7:
                add=add+1
            if numbersHyb[iCount][0]>6:
                exclusHybUp.append(exclusSho[iCount+add][0:35])
                for jCount in range(len(numbersHyb[iCount])):
                    empty=""
                    lenDigit=len(str(numbersHyb[iCount][jCount]))
                    # Determining how many spaces must be put in list
                    if jCount==0:
                        spaces = 4-lenDigit
                    else:
                        spaces = 5-lenDigit
                    for kCount in range(spaces):
                        empty=empty+" "
                    exclusHybUp[iCount] = (exclusHybUp[iCount] + empty
                        + str(numbersHyb[iCount][jCount]))
                exclusHybUp[iCount]=exclusHybUp[iCount]+"\n"
            else:
                exclusHybUp.append(exclusSho[iCount+add][0:35])
                for jCount in range(len(numbersHyb[iCount])):
                    empty=""
                    lenDigit=len(str(numbersHyb[iCount][jCount]))
                    # Determine how many spaces must be put in list
                    spaces = 5-lenDigit
                    if jCount==0:
                        spaces = 4-lenDigit
                    else:
                        spaces = 5-lenDigit
                    for kCount in range(spaces):
                        empty=empty+" "
                    exclusHybUp[iCount] = (exclusHybUp[iCount] + empty
                        + str(numbersHyb[iCount][jCount]))
                exclusHybUp[iCount]=exclusHybUp[iCount]+"\n"
    # Breaking all exceeding lines
    emptyLine=""
    # Spaces for new line to add
    for iCount in range(39):
        emptyLine=emptyLine+" "
    addNum = 0
    for jCount in range(len(numbersHyb)):
        multip=2
        start=7
        # Looping over how many times 6 fits into the numbers of exclusions,
        # i.e. how many new lines do we need.  
        for iCount in range(int(numbersHyb[jCount][0]/6.0)):
            if start==numbersHyb[jCount][0]+1:
                break
            exclusHybUp.insert(jCount+addNum+1,emptyLine)
            # From first line of the exclusions, just keep the first 6
            # exclusions.  
            exclusHybUp[jCount+addNum]=exclusHybUp[jCount+addNum][0:69]+"\n"
            for lCount in range(start,numbersHyb[jCount][0]+1):
                empty=""
                lenDigit=len(str(numbersHyb[jCount][lCount]))
                # Determining how many spaces must be put in list
                spaces = 5-lenDigit
                for kCount in range(spaces):
                    empty=empty+" "
                exclusHybUp[jCount+addNum+1] = (exclusHybUp[jCount+addNum+1]
                    + empty + str(numbersHyb[jCount][lCount]))
                if lCount>=6*multip:
                    start=start+6
                    multip=multip+1
                    break
            exclusHybUp[jCount+addNum+1] = exclusHybUp[jCount+addNum+1] + "\n"
            addNum=addNum+1
    return exclusHybUp

# Adds atom numbers of hybrid to exclusion list of newly inserted atoms.  
# lastNumLo is the atom number of the last atom in longer residue,
# analogous lastNumSho.  addAtom is the amount of atoms, that will be added
# to hybrid mtb.  
def updateAtom(atomInfo, lastNumSho, lastNumLo, hypPro):
    addAtom = lastNumSho-4
    addEx = range(5+addAtom,addAtom+lastNumLo+1)
    addEx.insert(0,len(addEx))
    numbersHyb=[]
    for iCount in range(addAtom):
        numbersHyb.append([])
        for jCount in range(len(addEx)):
            numbersHyb[iCount].append(addEx[jCount])
    return updateExclu(atomInfo, [], addAtom, False, numbersHyb, hypPro)

# Returns all bond lengths, angles, dihedrals... of new atoms.  In the first
# element list of returning list bondRes, there are saved 4 integers.  These
# are the updated #NB, #NBA, ...    
def updateBonds(atomInfo, bondsSho, bondsLo, withArg, withPro):
    bondRes=[[],[],[],[],[]]
    # Extracting atom numbers from atomInfo
    numbersAdd=[]
    if withArg:
        for iCount in range(len(atomInfo)):
            line = atomInfo[iCount].split()
            numbersAdd.append(int(line[0]))
    # Additional H atom must be considered when PRO/HYP is involved.  
    if withPro:
        numbersAdd.insert(0,2)
    # iCount iterates over lengths block, angle block...  
    found=False
    column=4
    for iCount in range(1,len(bondsSho)):
        # Check specifically the lengths (2 atoms define 1 length)
        if iCount==1:
            for jCount in range(len(numbersAdd)):
                # lCount runs over all lines in iCount-th block
                for lCount in range(len(bondsSho[iCount])):
                    # Fixating atom in atomInformation.  
                    # Check first and second entry of lCount-th line in bond
                    # length block of short residue.  
                    FELB=bondsSho[iCount][lCount][0]
                    SELB=bondsSho[iCount][lCount][1]
                    # If there is a defined length for an atom that must be
                    # added, adding it to bondRes.  
                    if numbersAdd[jCount]==FELB or numbersAdd[jCount]==SELB:
                        # Lengths are sorted by atom in first column.  
                        bondRes[iCount].append(bondsSho[iCount][lCount])
                        bondsLo[0][iCount-1]=bondsLo[0][iCount-1]+1
        if iCount==2:
            for jCount in range(len(numbersAdd)):
                # lCount runs over all lines in iCount-th block.  
                for lCount in range(len(bondsSho[iCount])):
                # Fixating atom in atomInformation.  
                # Check first, second and third entry of lCount-th line in
                # bond angles block of short residue.  
                    FELB=bondsSho[iCount][lCount][0]
                    SELB=bondsSho[iCount][lCount][1]
                    TELB=bondsSho[iCount][lCount][2]
                    # If there is a defined angle for an atom that must be
                    # added, adding it to bondRes.  
                    if (numbersAdd[jCount] == FELB
                            or numbersAdd[jCount]==SELB
                            or numbersAdd[jCount]==TELB):
                        bondRes[iCount].append(bondsSho[iCount][lCount])
                        bondsLo[0][iCount-1]=bondsLo[0][iCount-1]+1
        if iCount>2:
            for jCount in range(len(numbersAdd)):
                # lCount runs over all lines in iCount-th block.  
                for lCount in range(len(bondsSho[iCount])):
                    # Fixating atom in atomInformation. Checking first,
                    # second, third and fourth entry of lCount-th line in bond
                    # angles block of short residue.  
                    FELB=bondsSho[iCount][lCount][0]
                    SELB=bondsSho[iCount][lCount][1]
                    TELB=bondsSho[iCount][lCount][2]
                    FOELB=bondsSho[iCount][lCount][3]
                    # If there is a defined angle for an atom that must be
                    # added, adding it to bondRes.  
                    if (numbersAdd[jCount] == FELB 
                            or numbersAdd[jCount]==SELB 
                            or numbersAdd[jCount]==TELB 
                            or numbersAdd[jCount]==FOELB):
                        bondRes[iCount].append(bondsSho[iCount][lCount])
                        bondsLo[0][iCount-1]=bondsLo[0][iCount-1]+1
        # Eliminating duplicates in bond lengths
        for jCount in range(len(bondRes[iCount])):
            addBonds=0
            for kCount in range(jCount+1,len(bondRes[iCount])):
                if bondRes[iCount][jCount]==bondRes[iCount][kCount-addBonds]:
                    del bondRes[iCount][kCount-addBonds]
                    bondsLo[0][iCount-1]=bondsLo[0][iCount-1]-1
                    addBonds=addBonds+1
    for iCount in range(4):
        bondRes[0].append(bondsLo[0][iCount])
    return bondRes

def bondsToString(newBondInt,bondsStr,residues,shorterRes,withArg):
    newBonds=[]
    newBonds.append("# bonds\n")
    newBonds.append("#  NB\n")
    empty=" "
    spaces = len(str(newBondInt[0][1]))
    for kCount in range(4-spaces):
        empty=empty+" "
    newBonds.append(empty+str(newBondInt[0][0])+"\n")
    newBonds.append("#  IB   JB  MCB\n")
    for iCount in range(len(bondsStr[0])):
        newBonds.append(bondsStr[0][iCount])
    newBonds.append("### "+residues[shorterRes][1]+" bonds\n")
    for iCount in range(len(newBondInt[1])):
        for jCount in range(len(newBondInt[1][iCount])):
            empty=" "
            spaces = len(str(newBondInt[1][iCount][jCount]))
            if jCount==0:
                for kCount in range(4-spaces):
                    empty=empty+" "
            else:
                for kCount in range(4-spaces):
                    empty=empty+" "
            if jCount==len(newBondInt[1][iCount])-1:
                newBonds.append(empty+str(newBondInt[1][iCount][jCount])+"\n")
            else:
                newBonds.append(empty+str(newBondInt[1][iCount][jCount]))
    newBonds.append("# bond angles\n")
    newBonds.append("# NBA\n")
    empty=" "
    spaces = len(str(newBondInt[0][1]))
    for kCount in range(4-spaces):
        empty=empty+" "
    newBonds.append(empty+str(newBondInt[0][1])+"\n")
    newBonds.append("#  IB   JB   KB  MCB\n")
    for iCount in range(len(bondsStr[1])):
        newBonds.append(bondsStr[1][iCount])
    newBonds.append("### "+residues[shorterRes][1]+" bond angles\n")
    for iCount in range(len(newBondInt[2])):
        for jCount in range(len(newBondInt[2][iCount])):
            empty=" "
            spaces = len(str(newBondInt[2][iCount][jCount]))
            if jCount==0:
                for kCount in range(4-spaces):
                    empty=empty+" "
            else:
                for kCount in range(4-spaces):
                    empty=empty+" "
            if jCount==len(newBondInt[2][iCount])-1:
                newBonds.append(empty+str(newBondInt[2][iCount][jCount])+"\n")
            else:
                newBonds.append(empty+str(newBondInt[2][iCount][jCount]))
    newBonds.append("# improper dihedrals\n")
    newBonds.append("# NIDA\n")
    empty=" "
    spaces = len(str(newBondInt[0][2]))
    for kCount in range(4-spaces):
        empty=empty+" "
    newBonds.append(empty+str(newBondInt[0][2])+"\n")
    newBonds.append("#  IB   JB   KB   LB  MCB\n")
    for iCount in range(len(bondsStr[2])):
        newBonds.append(bondsStr[2][iCount])
    newBonds.append("### "+residues[shorterRes][1]+" improper dihedrals\n")
    for iCount in range(len(newBondInt[3])):
        for jCount in range(len(newBondInt[3][iCount])):
            empty=" "
            spaces = len(str(newBondInt[3][iCount][jCount]))
            if jCount==0:
                for kCount in range(4-spaces):
                    empty=empty+" "
            else:
                for kCount in range(4-spaces):
                    empty=empty+" "
            if jCount==len(newBondInt[3][iCount])-1:
                newBonds.append(empty+str(newBondInt[3][iCount][jCount])+"\n")
            else:
                newBonds.append(empty+str(newBondInt[3][iCount][jCount]))
    newBonds.append("# dihedrals\n")
    newBonds.append("# NDA\n")
    empty=" "
    spaces = len(str(newBondInt[0][3]))
    for kCount in range(4-spaces):
        empty=empty+" "
    newBonds.append(empty+str(newBondInt[0][3])+"\n")
    newBonds.append("#  IB   JB   KB   LB  MCB\n")
    for iCount in range(len(bondsStr[3])):
        newBonds.append(bondsStr[3][iCount])
    newBonds.append("### "+residues[shorterRes][1]+" dihedrals\n")
    for iCount in range(len(newBondInt[4])):
        for jCount in range(len(newBondInt[4][iCount])):
            empty=" "
            spaces = len(str(newBondInt[4][iCount][jCount]))
            for kCount in range(4-spaces):
                empty=empty+" "
            if jCount==len(newBondInt[4][iCount])-1:
                newBonds.append(empty+str(newBondInt[4][iCount][jCount])+"\n")
            else:
                newBonds.append(empty+str(newBondInt[4][iCount][jCount]))
    return newBonds

# Writes the modified.mtb file with all atoms.  
def atomsToFile(exclusHybUp, excludedAtoms, newBonds,
                contentsMtb, residues, shorterRes,
                longerRes, outFileName, lastNumSho,
                replace, withArg, hypPro):
    eraseFrom = -1
    eraseTo = -1
    oldLength = len(contentsMtb)
    newLength = (len(contentsMtb) + len(exclusHybUp)
                 + len(excludedAtoms) + len(newBonds))
    sizeRep = len(replace)
    if hypPro:
        searchC = "    1 N "
    else:
        searchC = "    3 CA"
    iCount = 0
    while hypPro and "#ATOM  " != contentsMtb[iCount][0:7]:
        iCount = iCount + 1
    for iCount in range(newLength):
        if eraseTo<0 and "# RNME\n" == contentsMtb[iCount]:
            contentsMtb[iCount+1] = (residues[shorterRes][0] + "_"
            + residues[longerRes][0]+"\n")
        if searchC == contentsMtb[iCount][0:8]:
            # eraseFrom points to 3 CA line
            eraseFrom = iCount
            for jCount in range(newLength):
                if "    5 New"==contentsMtb[jCount][0:9]:
                    # eraseTo points to 5 New line
                    eraseTo=jCount
                    break
                elif jCount==oldLength-1:
                    ### CAREFUL: eraseTo = eraseFrom + 2 fixed ALA
                    if withArg:
                        eraseTo=eraseFrom + 2
                    break
            for jCount in range(eraseFrom,eraseTo):
                del contentsMtb[eraseFrom]
            if withArg:
                for jCount in range(len(exclusHybUp)):
                    contentsMtb.insert(eraseFrom + jCount,exclusHybUp[jCount])
            # eraseFrom points to 5 New line
            if not hypPro:
                eraseFrom=iCount+len(exclusHybUp)
            else:
                eraseFrom=iCount+7
            contentsMtb.insert(eraseFrom,("### "+residues[shorterRes][1]
                                          +" atoms\n"))
            while "New" in contentsMtb[eraseFrom+1]:
                del contentsMtb[eraseFrom+1]
            for jCount in range(len(excludedAtoms)):
                contentsMtb.insert(eraseFrom+jCount+1,excludedAtoms[jCount])
            # eraseFrom points to first atom, that already was in hybrid
            if not hypPro:
                eraseFrom=iCount+len(exclusHybUp)+len(excludedAtoms)+2
            else:
                eraseFrom=iCount+len(exclusHybUp)+len(excludedAtoms)+9
            contentsMtb.insert(eraseFrom-1,"###\n")
            # Adding an 'A' to all atom name clashes of hybrid.  
            # kCount runs over contentsMtb.  
            kCount=eraseFrom
            while(sizeRep>0):
                line=contentsMtb[kCount].split()
                if len(line)>1 and line[1]==replace[0][1][0:-1]:
                    replace[0][1]=replace[0][1][0:-1]+"A"
                    spaces=len(replace[0][1])
                    empty=" "
                    for lCount in range(6-spaces):
                        empty=empty+" "
                    contentsMtb[kCount] = (contentsMtb[kCount][0:6] 
                          +replace[0][1] + empty + contentsMtb[kCount][13:])
                    del replace[0]
                    sizeRep=sizeRep-1
                    kCount=eraseFrom
                else:
                    kCount=kCount+1
            for jCount in range(newLength):
                if "# bonds\n"==contentsMtb[jCount]:
                    # eraseFrom points to '# bonds' line
                    eraseFrom=jCount
                    break
            for jCount in range(newLength):
                if "# LJ exceptions\n"==contentsMtb[jCount]:
                    # eraseTo points to #LJ exceptions line
                    eraseTo=jCount
                    break
            for jCount in range(eraseTo-eraseFrom):
                del contentsMtb[eraseFrom]
            for jCount in range(len(newBonds)):
                contentsMtb.insert(eraseFrom+jCount,newBonds[jCount])
            break
    outFile = open((residues[shorterRes][1] + "_" + residues[longerRes][1]
            + "_modified.mtb"),'w')
    for iCount in range(len(contentsMtb)):
        outFile.write(contentsMtb[iCount])
    outFile.close()

# Decides, whether CB should be in .ptp file and writes one .ptp file for each direction.  
def perturb(
                massFileName, residues, forceLines,
                shorterRes, longerRes, lastNumSho,
                lastNumLo, exclusSho, exclusHyb,
                rememberAI, atomInfoHyb, renamedAtoms,
                withArg, withPro, hypPro,
                withAla, withGly):
    massFile = open(massFileName, 'r')
    # Reading in the .mtb file to be changed and closing the file
    masses = massFile.readlines()
    massFile.close()
    # Needed for reading lines of modified .mtb file
    if withPro:
        addName = "_2"
    else:
        addName = ""
    readModif = (residues[shorterRes][1] + "_"
                + residues[longerRes][1] + "_modified" + addName + ".mtb")
    ptp1 = open((residues[shorterRes][1] + "_"
                + residues[longerRes][1] + ".ptp"), 'w')
    ptp2 = open((residues[longerRes][1] + "_"
                + residues[shorterRes][1] + ".ptp"), 'w')
    ptp1name = (residues[shorterRes][1] + "_"
                + residues[longerRes][1] + ".ptp")
    ptp2name = (residues[longerRes][1] + "_"
                + residues[shorterRes][1] + ".ptp") 
    title1 = ("template perturbation topology for perturbation of "
              + residues[shorterRes][2].lower() + " (state A) to "
              + residues[longerRes][2].lower() + " (state B)\n")
    title2 = ("template perturbation topology for perturbation of "
              + residues[longerRes][2].lower() + " (state A) to "
              + residues[shorterRes][2].lower() + " (state B)\n")

    # Checking if CB must be part of perturbation
    lineSho = exclusSho[1].split()
    lineLo = exclusHyb[1].split()
    lineSho = lineSho[1:5]
    lineLo = lineLo[1:5]
    equal = True
    case = 0
    iCount = 0
    while equal and iCount < len(lineSho):
        if lineSho[iCount]!=lineLo[iCount]:
            equal = False
        iCount = iCount+1
    if hypPro:
        # PRO <-> HYP
        equal = True
        subtract = 5
        case = 7
    else:
        if withPro:
            if not (withGly or withAla):
                # X \ {GLY, ALA} <-> PRO, HYB
                subtract = -1
                equal = False
                case = 4
            elif withGly:
                # GLY <-> PRO OR HYP
                subtract = -1
                equal = False
                case = 5
            else:                
                # ALA <-> PRO OR HYP
                subtract = -1
                equal = False
                case = 6
        else:
            if withGly:
                if withAla:
                    # ALA <-> GLY
                    subtract = 2
                    case = 8
                else:
                    # GLY <-> X \ {PRO, HYP, ALA}
                    subtract = 2
                    equal = False
                    case = 2
            else:
                if withAla:
                    # ALA <-> X \ {PRO, HYP, GLY}
                    subtract = 4
                    equal = False
                    case = 3
                else:
                    # X \ {PRO, HYP, GLY, ALA} <-> X \ {PRO, HYP, GLY, ALA}
                    subtract = 4
                    case = 1
    addAtom = lastNumSho-subtract
    if equal:
        addEx = range(1,lastNumLo+addAtom-3)
    else:
        addEx = range(1,lastNumLo+addAtom-2)
    res = "X"
    text1 = "TITLE\n"
    text2 = "END\nPERTATOMPARAM\n# number of perturbed atoms\n"
    # Number of atom to perturb
    text3 = getEmpty(4,len(str(len(addEx))))+str(len(addEx))
    text4 = ("\n#\n# NR RES NAME IAC(A) MASS(A)  "
             + "CHARGE(A) IAC(B) MASS(B)  CHARGE(B)   ALJ  ACRF\n")
    # Determining the actual perturbation atoms and gets the whole string
    # to insert in first .ptp file.  
    firstDirection = True
    text5 = getPertList(masses, readModif, addEx,
                        addAtom, exclusSho, exclusHyb,
                        renamedAtoms, equal, firstDirection,
                        withArg, withPro, hypPro, withAla, withGly, case)
    firstDirection = False
    text6 = getPertList(masses, readModif, addEx, 
                        addAtom, exclusSho, exclusHyb,
                        renamedAtoms, equal, firstDirection,
                        withArg, withPro, hypPro, withAla, withGly, case)
    text7 = "END"
    ptp1.write(text1)
    ptp1.write(title1)
    ptp1.write(text2)
    ptp1.write(text3)
    ptp1.write(text4)
    for iCount in text5:
        ptp1.write(iCount),
    ptp1.write(text7)

    ptp2.write(text1)
    ptp2.write(title2)
    ptp2.write(text2)
    ptp2.write(text3)
    ptp2.write(text4)
    for iCount in text6:
        ptp2.write(iCount),
    ptp2.write(text7)

    ptp1.close()
    ptp2.close()
    return case

# Determines the actual perturbation atoms by opening modified.mtb and gets
# the whole string to insert in first or second .ptp file.  exclusX holds
# information of CA and CB.  
def getPertList(masses, readModif, addEx, 
                addAtom, exclusSho, exclusHyb,
                renamedAtoms, equal, firstDirection,
                withArg, withPro, hypPro, withAla, withGly, case):
    sidechainF = open(readModif,'r')
    sidechain = sidechainF.readlines()
    sidechainF.close()
    
    # Deleting all lines of modified.mtb up to first side chain atom line
    jCount = 0
    found = False
    if case == 1 or case == 3:
        search = "###"
        searchInd = 0
    elif case == 2 or case == 8:
        search = "CA"
        searchInd = 1
    elif case == 7:
        search = "CG"
        searchInd = 1
    elif case == 4 or case == 5 or case == 6:
        search = "N"
        searchInd = 1
    while jCount <= len(sidechain) and not found:	
        if (len(sidechain[jCount].split())<3 
            or sidechain[jCount].split()[searchInd] != search):
            jCount = jCount + 1
        else:
            found = True
    for iCount in range(jCount+1):
        del sidechain[0]
    found = False
    jCount = 0
    # Deleting all lines of modified.mtb from first line after last side 
    # chain atom.  
    while jCount <= len(sidechain) and not found:
        if (len(sidechain[jCount].split()) == 1 
               and sidechain[jCount].split()[0] == "###"):
            del sidechain[jCount]
        elif (len(sidechain[jCount].split()) < 3
               or (sidechain[jCount].split()[0] != "###" 
               and sidechain[jCount].split()[2] != "atoms")):
            jCount = jCount + 1
        else:
            found = True
    if hypPro or withArg and not withPro:
        for iCount in range(len(sidechain)-jCount):
            del sidechain[jCount]
    jCount = 0
    lenSide = len(sidechain)
    found = False
    while withPro and jCount < lenSide and not found:
        checkLine = sidechain[jCount].split()
        if len(checkLine) > 1:
            if checkLine[1] == "trailing":
                found = True
            else:
                jCount = jCount + 1
        else:
            jCount = jCount + 1
    if withPro:
        for iCount in range(len(sidechain)-jCount):
            del sidechain[jCount]
    jCount = 0
    lenSide = len(sidechain)
    # Deleting all lines, that have less than 7 entries, so that first line
    # stores first atom, second line stores second atom and so on.  
    # All other lines are irrelevant.  
    while jCount < lenSide:
        if len(sidechain[jCount].split())<7:
            del sidechain[jCount]
            lenSide = lenSide - 1
        else:
            jCount = jCount + 1
    pertList = []
    exclus = []
    exclus.append(exclusSho)
    exclus.append(exclusHyb)
    indent = 1
    jCount = 0
    while repr("# N     ATMAS   ATMASN\n") != repr(masses[jCount]):
        jCount = jCount + 1
    # jCount+1 points on mass of hydrogen
    hydrogen = jCount + 1
    # Adding CB atom as first atom, if they're different in both residues
    if not equal and not withPro:
        nextLine=""
        half=""
        temp=""
        for iCount in range(2):
            alphaBeta = copy.deepcopy(exclus[iCount])
            # Atom masses of CB/CA atoms.  
            if case == 1 or case == 3 or case == 4 or case == 6:
                subtract = 0
                backC = "CB"
            elif case == 2 or case == 5 or case == 8:
                subtract = 1
                backC = "CA"

            massCB = alphaBeta[1-subtract].split()[3]
            jCount = hydrogen
            # All masses are stored up to line 70 in .ifp file.  
            # Searching for right line.  
            while (jCount <= 70 
                    and int(massCB) != int(masses[jCount].split()[0])):
                jCount = jCount+1
            # exclusX[1].split()[1] is the renamed atom name,
            # ...[2] is the IAC, ...[4] is the charge. Building the line.  
            if iCount==0:
                prefix = (getEmpty(3,1) + str(1) + "   X  " + backC
                       + getEmpty(5,len(backC)))
                half = (" " + alphaBeta[1-subtract].split()[2] + "   "
                     + getEmpty(7, len(masses[jCount].split()[1]))
                     + masses[jCount].split()[1] + "   "
                     + getEmpty(8,len(alphaBeta[1-subtract].split()[4]))
                     + alphaBeta[1-subtract].split()[4])
            else:
                half = (" " + alphaBeta[1-subtract].split()[2] + "   "
                     + getEmpty(7, len(masses[jCount].split()[1]))
                     + masses[jCount].split()[1] + "   "
                     + getEmpty(8, len(alphaBeta[1-subtract].split()[4]))
                     + alphaBeta[1-subtract].split()[4])
            if firstDirection:
                if iCount==0:
                    nextLine = prefix + nextLine + half
                else:
                    nextLine = nextLine + "  "+half
            else:
                if iCount==0:
                    temp = copy.deepcopy(half)
                else:
                    nextLine = prefix + half + "  " + temp
        nextLine = (nextLine + getEmpty(7, len("1.0")) + "1.0"
                 + getEmpty(5, len("1.0")) + "1.0\n")
        pertList.append(nextLine)
    elif withPro and firstDirection:
        pertList.insert(0,("  1   X  N      6   14.0067   -0.31000    6"
                           + "   14.0067    0.00000    1.0  1.0\n"))
    elif withPro and not firstDirection:
        pertList.insert(0,("  1   X  N      6   14.0067    0.00000    6"
                           + "   14.0067   -0.31000    1.0  1.0\n"))
    elif hypPro and firstDirection:
        pertList.insert(0,("  1   X  CG    18    14.027    0.00000   14"
                           + "    13.019    0.26600    1.0  1.0\n"))
    elif hypPro and not firstDirection:
        pertList.insert(0,("  1   X  CG    14    13.019    0.26600   18"
                           + "    14.027    0.00000    1.0  1.0\n"))
    else:
        indent = 0
    if case == 2 or case == 7:
        indent = 1
        addAtom = 0
    elif case == 3:
        indent = 1
    elif case == 8:
        addAtom = 0
    elif case == 6 or case == 5 or case == 4:
        indent = 1
        addAtom = 1
    switched = False
    iCount = 0
    while iCount < len(sidechain):
        if "#" == sidechain[iCount].split()[0]:
            del sidechain[iCount]
        iCount = iCount + 1
    for iCount in range(len(addEx)-indent):
        if withPro and not switched and iCount >= len(exclus[1]):
            firstDirection = not firstDirection
            switched = True
        nextLine=""
        sideChainLine = copy.deepcopy(sidechain[iCount])
        # Atom masses of atoms. Remember: there are no changes in masses in 
        # 'regular side chain atoms' (i.e. non-CB atoms).  
        massNum = sideChainLine.split()[3]
        jCount = hydrogen
        # All masses are stored up to line 70 in .ifp file.  
        # Searching for right line.
        while jCount <= 70 and int(massNum) != int(masses[jCount].split()[0]):
            jCount = jCount+1
        mass = masses[jCount].split()[1]
        lj = "1.0"
        # kCount runs over first and second half of each 
        # line (active/inactive state).  
        for kCount in range(2):
            if kCount==0:
                num = addEx[iCount+indent]
                name = sideChainLine.split()[1]
                if iCount < addAtom and firstDirection:
                    iac = sideChainLine.split()[2]
                    charge = sideChainLine.split()[4]
                elif iCount < addAtom and not firstDirection:
                    iac = "22"
                    charge = "0.00000"
                elif iCount >= addAtom and firstDirection:
                    iac = "22"
                    charge = "0.00000"
                else:
                    iac = sideChainLine.split()[2]
                    charge = sideChainLine.split()[4]
                half = (getEmpty(3, len(str(num))) + str(num) + "   X  "
                     + name + getEmpty(5, len(name)) + " "
                     + getEmpty(2, len(iac)) + iac + "   "
                     + getEmpty(7, len(mass)) + mass + "   "
                     + getEmpty(8,len(charge)) + charge)
            else:
                if iCount < addAtom and firstDirection:
                    iac = "22"
                    charge = "0.00000"
                elif iCount < addAtom and not firstDirection:
                    iac = sideChainLine.split()[2]
                    charge = sideChainLine.split()[4]
                elif iCount >= addAtom and firstDirection:
                    iac = sideChainLine.split()[2]
                    charge = sideChainLine.split()[4]
                else:
                    iac = "22"
                    charge = "0.00000"
                half = ("   " + getEmpty(2, len(iac)) + iac + "   "
                     + getEmpty(7, len(mass)) + mass + "   "
                     + getEmpty(8,len(charge)) + charge)
            nextLine = nextLine + half
        nextLine = (nextLine + getEmpty(7, len(lj))+ lj + getEmpty(5,len(lj)) 
                 + lj + "\n")
        pertList.append(nextLine)
    if withPro and firstDirection:
        del pertList[2]
        if case == 5:
            pertList.insert(2,("  3   X  CA    15    14.027    0.00000"
                            + "   14    13.019    0.00000    1.0  1.0\n"))
        else:
            pertList.insert(2,("  3   X  CA    14    13.019    0.00000"
                            + "   14    13.019    0.00000    1.0  1.0\n"))
    elif withPro and not firstDirection:
        del pertList[2]
        if case == 5:
            pertList.insert(2,("  3   X  CA    14    13.019    0.00000"
                            + "   15    14.027    0.00000    1.0  1.0\n"))
        else:
            pertList.insert(2,("  3   X  CA    14    13.019    0.00000"
                            + "   14    13.019    0.00000    1.0  1.0\n"))
    return pertList

# Returns an empty string, that depends on number of digits and length of
# empty prefix.  
def getEmpty(digits, empDig):
    empty=""
    for iCount in range(digits-empDig):
        empty = empty + " "
    return empty

def header(shorterRes, longerRes, residues, case):
    mod2 = [4, 5, 6]
    add = "" 
    if case in mod2:
        add = "_2"
    path = (residues[shorterRes][1] + "_" + residues[longerRes][1]
           + "_modified" + add + ".mtb")
    add1 = "# This is a hybrid building block to be used for perturbing\n"
    add2 = ("# " + residues[shorterRes][1] + " to "
            + residues[longerRes][1] + "\n")
    add3 = "# based on the GROMOS 54a7 force field.\n"
    addString = [add3, add2, add1]
    linesFile = open(path,'r+')
    lines = linesFile.readlines()
    while lines[0].split()[0] == "#":
        del lines[0]
    for iCount in addString:
        lines.insert(0,iCount)
    linesFile.seek(0)
    for iCount in lines:
        linesFile.write(iCount)
    linesFile.truncate()
    linesFile.close()

def rename(shorterRes, longerRes, residues, case):
    mod2 = [4, 5, 6]
    add = "" 
    if case in mod2:
        add = "_2"
    path = (residues[shorterRes][1] + "_" + residues[longerRes][1]
           + "_modified" + add + ".mtb")
    newName = (residues[shorterRes][1] + "_" + residues[longerRes][1]
           + ".mtb")
    os.rename(path,newName)

def delete(shorterRes, longerRes, residues):
    delete = [(residues[shorterRes][1]
                  + "_" + residues[longerRes][1] + "_initial.mtb"),
              (residues[shorterRes][1]
                  + "_" + residues[longerRes][1] + "_initial_1.mtb"),
              (residues[shorterRes][1]
                  + "_" + residues[longerRes][1] + "_initial_2.mtb"),
              (residues[shorterRes][1]
                  + "_" + residues[longerRes][1] + "_modified_1.mtb"),
              "add_atom_1.arg", "add_atom_2.arg", "add_atom.arg"]
    for iCount in delete:
        try:
            os.remove(iCount)
        except OSError:
            pass

def printAll(atomInfo, one, rememberAI,
             two, exclusSho, three,
             atomInfoHyb, four, exclusHyb,
             five, exclusHybUp, six,
             atomInfoEx, seven, bonds,
             eight, column1, oldAmount,
             newBonds, nine):
    if one:
        print "Atoms of short residue after renaming:"
        for iCount in range(len(atomInfo)):
            print atomInfo[iCount],
    if two:
        print "Atoms of short residue before renaming:"
        for iCount in range(len(rememberAI)):
            print rememberAI[iCount],
    if three:
        print "Exclusions short residue:"
        for iCount in range(len(exclusSho)):
            print exclusSho[iCount],
    if four:
        print "Atoms of hybrid before updating:"
        for iCount in range(len(atomInfoHyb)):
            print atomInfoHyb[iCount],
    if five:
        print "Exclusions hybrid:"
        for iCount in range(len(exclusHyb)):
            print exclusHyb[iCount],
    if six:
        print "Numbers update exclusion list of hybrid:"
        for iCount in range(len(exclusHybUp)):
            print exclusHybUp[iCount],
    if seven:
        print "Updated exclusion list with hybrid chain atoms:"
        for iCount in range(len(atomInfoEx)):
            print atomInfoEx[iCount],
    if eight:
        for iCount in range(len(bonds[column1])):
            if iCount==0:
                if column1==1:
                    print (str(bonds[0][column1-1] - oldAmount[column1-1])
                        +" hybrid bond lengths added:")
                elif column1==2:
                    print (str(bonds[0][column1-1] - oldAmount[column1-1])
                        +" hybrid bond angles added:")
                elif column1==3:
                    print (str(bonds[0][column1-1] - oldAmount[column1-1])
                    + " hybrid improper dihedrals added:")
                elif column1==4:
                    print (str(bonds[0][column1-1] - oldAmount[column1-1])
                    + " hybrid dihedrals before added:")
            print bonds[column1][iCount]
    if nine:
        for iCount in range(len(newBonds)):
            print newBonds[iCount],

def printFiles(shorterRes, longerRes, residues, one,
               head, two, three, case, 
               four):
    mod2 = [4, 5, 6]
    add = "" 
    if case in mod2:
        add = "_2"
    ptp1 = residues[shorterRes][1] + "_" + residues[longerRes][1] + ".ptp"
    ptp2 = residues[longerRes][1] + "_" + residues[shorterRes][1] + ".ptp"
    path = (residues[shorterRes][1] + "_" + residues[longerRes][1]
           + "_modified" + add + ".mtb")
    if one:
        print "\n"+path
        mtbFile = open(path,'r')
        mtbLines = mtbFile.readlines()
        mtbFile.close()
        print "Topology:"
        for iCount in range(head):
            print "\t"+mtbLines[iCount],
        print
    if two:
        print "\n"+ptp1
        ptp1File = open(ptp1,'r')
        ptp1Lines = ptp1File.readlines()
        ptp1File.close()
        print "Perturbation forward:"
        for iCount in ptp1Lines:
            print "\t"+iCount,
        print
    if three:
        print "\n"+ptp2
        ptp2File = open(ptp2,'r')
        ptp2Lines = ptp2File.readlines()
        ptp2File.close()
        print "Perturbation backwards:"
        for iCount in ptp2Lines:
            print "\t"+iCount,
        print
    if four:
        print "Case: " + str(case)

main(sys.argv)
