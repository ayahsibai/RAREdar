"""
Created on Wed May 9 12:08:07 2022

@author: Frank Zhuang
"""

import utility as u
import re as r


def main():

    Motif = r'ACT[AGT]G[AGT]' #This is the original, 5'-3' 8 variants of RARE motif
    Complement = r'TGA[ACT]C[ACT]' #This is the COMPLEMENT variants to the RARE motif (3'-5') on the ANTISENSE strand
    GeneDictionary = u.readFasta('Retinal_Candidate_Set/D_rerio_Retinal_Sequence.txt')
    #print(GeneDictionary)

    hitDictionary,positionDictionary,sequenceDictionary = RAREdar(GeneDictionary, Motif, Complement) #Run RAREdar

    truepositionDictionary = con_coord(positionDictionary) #RAREdar outputs coordinates RELATIVE to the gene of interest. This function converts the data to TRUE coordinates on the chromosome

    auto_merger(Motif, Complement, hitDictionary, truepositionDictionary, sequenceDictionary) #Outputs raw data as a tab delimited file

    return

'''
#RAREdar
#Read along sequences and compare every 6-bp window to every 6-mer motif in the
list of motifs, reports all instances when a 6-bp window matches a predescribed
motif and have a DR 5bps after it
#Input: Dictionary with gene reference number as keys and gene sequence as values
#Input: A list of all 6-mer motifs to be compared against
#Output: Dictionary with gene reference number as keys and number of hits, the
position of every hit, and the sequence of every hit as values.
'''

def RAREdar(dictionary, DRList, RDRList):

    hits = 0 #define object and assign initial value for hits, positions and sequence
    hitDictionary = {}
    positionDictionary = {}
    hitPosition = []
    sequenceDictionary = {}
    hitSequence = []

    DRReverse = r'[TGA]G[TGA]TCA' #Create list of reversed version of RAREs
    RDRReverse = r'[TCA]C[TCA]AGT'

    for name in dictionary.keys():
        title = name.split()
        sequence = dictionary[name] #assign the sequence of a gene as a string object
        for bp in range(len(sequence)-17):
            window = sequence[bp:bp+6] #create 6-bp window that reads along the sequence string
            repeat = sequence[bp+11:bp+17]
            if r.search(DRList, window) and window == repeat: #compare every window to every 6-mer in CODING RARE list
                hits = hits + 1 #record valid hit
                hitPosition = hitPosition + [bp]
                hitSequence = hitSequence + [sequence[bp:bp+17]]
            elif r.search(RDRList, window) and window == repeat: #compare every window to every 6-mer in COMPLEMENTARY RARE list
                hits = hits + 1 #record valid hit
                hitPosition = hitPosition + [bp]
                hitSequence = hitSequence + [sequence[bp:bp+17]]
            elif r.search(DRReverse, window) and window == repeat: #compare every window to every 6-mer in REVERSED RARE list
                hits = hits + 1 #record valid hit
                hitPosition = hitPosition + [bp]
                hitSequence = hitSequence + [sequence[bp:bp+17]]
            elif r.search(RDRReverse, window)and window == repeat: #compare every window to every 6-mer in REVERSED COMPLEMENTARY RARE list
                hits = hits + 1 #record valid hit
                hitPosition = hitPosition + [bp]
                hitSequence = hitSequence + [sequence[bp:bp+17]]
        if hits != 0:
            hitDictionary[name] = hits #Assign hit info as values to the gene reference number as keys
            positionDictionary[name] = hitPosition
            sequenceDictionary[name] = hitSequence
        hits = 0 #reinitialize all objects
        hitPosition = []
        hitSequence = []

    return hitDictionary, positionDictionary, sequenceDictionary

'''
#CONCOORD
#Takes relative space coordinates of found RAREs and convert them to true coordinates
#Input :: a dictionary of gene name to relative coordinates
#Output :: a dictionary of gene name to true coordinates
'''

def con_coord(dictionary):
    newdict = {}
    for name in dictionary.keys():
        title = name.split(' ') #Add ID info to every RARE hit
        coordRange = title[1][6:]
        position = coordRange.split(':')
        chromosome = position[0]
        coordinate = position[1].split('-')
        start = int(coordinate[0])
        end = int(coordinate[1])
        status = title[4][-1]
        truecoord = []
        if status == '+': #Situation for gene on forward strand
            for coord in dictionary[name]:
                truecoord = truecoord + [start+coord] #True coordinate by adding relative coordinate to gene starting coordinate
        elif status == '-': #Situation for gene on reversed starnd
            for coord in dictionary[name]:
                truecoord = truecoord + [end-coord-17] #True coordinate by subtracting relative coordinate to gene ending coordinate
        newdict[name] = truecoord
    return newdict

'''
#AUTO MERGER
#Takes all output from RAREdar and create a single, tab-delimited file of all data types
#Input :: 3 dictionaries, representing hit count, position and sequence
#Output :: tab-delimited plain text file with position and sequence data
'''

def auto_merger(DRList, RDRList, hitDictionary, positionDictionary, sequenceDictionary):
    DRReverse = r'[TGA]G[TGA]TCA' 
    RDRReverse = r'[TCA]C[TCA]AGT'

    fileName = 'output/RAREdar_Results_Retinal.txt'
    fileWrite = open(fileName,'w')
    fileWrite.write('Chromosome'+'\t'+'Gene'+'\t'+'Mode'+'\t'+'Coordinate'+'\t'+'Original Sequence'+'\t'+'Sense Sequence'+'\n')
    Entry = ''

    seenCoordinates = set()  # Set to store seen coordinates

    for name in hitDictionary.keys():
        title = name.split(' ')
        coordRange = title[1][6:]
        position = coordRange.split(':')
        chromosome = position[0]
        geneName = title[0]
        coordList = positionDictionary[name]
        sequenceList = sequenceDictionary[name]

        for i in range(len(coordList)):
            coordinate = coordList[i]
            sequence = sequenceList[i]

            # Check if the coordinate has been seen before
            if coordinate in seenCoordinates:
                continue  # Skip the entry if the coordinate is not unique

            seenCoordinates.add(coordinate)  # Add the coordinate to seenCoordinates set

            Entry = ''
            if r.search(DRList, sequence[0:6]):
                mode = 'Forward Coding'
                senseSequence = str(sequence)
            elif r.search(RDRList, sequence[0:6]):
                mode = 'Forward Complement'
                senseSequence = ''.join(['T' if bp == 'A' else 'A' if bp == 'T' else 'C' if bp == 'G' else 'G' for bp in sequence])
            elif r.search(DRReverse, sequence[0:6]):
                mode = 'Reversed Coding'
                senseSequence = str(sequence[::-1])
            elif r.search(RDRReverse, sequence[0:6]):
                mode = 'Reversed Complement'
                revSequence = str(sequence[::-1])
                senseSequence = ''.join(['T' if bp == 'A' else 'A' if bp == 'T' else 'C' if bp == 'G' else 'G' for bp in revSequence])
            else:
                mode = 'Exception'
                senseSequence = 'n/a'
            
            Entry = Entry + chromosome +'\t' + geneName +'\t' + mode +'\t' + str(coordinate) +'\t' + sequence + '\t' + senseSequence + '\n'
            fileWrite.write(Entry)
    return

'''
#DR SLIDER
#Use a target pattern and a gap width to find all matching hit of this target within a sequence
#Input :: a target pattern, a genetic sequence, the length of gap between direct repeats
#Output :: number of hits, list of position of all hits, list of sequence of all hits
'''

def dr_slider(target, sequence, gap):
    hits = 0
    hitSequence = []
    hitPosition = []
    targetLength = len(target)*2+gap
    for bp in range(len(sequence)-targetLength):
        window = sequence[bp:bp+len(target)] #create window that reads along the sequence string
        repeat = sequence[bp+targetLength-len(target):bp+targetLength]
        if window == repeat and window == target:
            hits = hits + 1 #record valid hit
            hitPosition = hitPosition + [bp]
            hitSequence = hitSequence + [sequence[bp:bp+17]]
    return hits, hitPosition, hitSequence

main()

