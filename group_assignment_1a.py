# Group Assignment 1a - Team 5
# Tyler Brian, Gracie Lovell, Ella Stone

# This function applies the input validation that makes sure that the given sequence is a valid DNA sequence
def InputValidation(DNA: str) -> bool:
    for c in DNA:
        if c not in ("A", "C", "G", "T"):
            return False
    return True

# This function finds where the consensus sequences are in the provided DNA sequence
def FindIndexOfCSeqs(DNA: str):
    start_motif = "TTGACA"  # start of promoter (-35 CS)
    end_motif = "TATAAT"  # end of promoter (-10 CS)
    indexOfPotentialStartPromoters = []
    indexOfPotentialEndOfPromoters = []
    dnaLength = len(DNA)
    for i in range(0, dnaLength - 5):
        segment = DNA[i:i+6]
        if segment == start_motif:
            indexOfPotentialStartPromoters.append(i)
        if segment == end_motif:
            indexOfPotentialEndOfPromoters.append(i)
    return indexOfPotentialStartPromoters, indexOfPotentialEndOfPromoters

# This function makes sure that the promoter is valid
# It is valid if the end sequence comes after the start sequence and there are 16–19 nucleotides between them
def ValidateCSeqs(starts, ends):
    validPromoters = []
    for s in starts:
        for e in ends:
            spacer = e - (s + 6)
            if 16 <= spacer <= 19:
                validPromoters.append((s, e))
    return validPromoters

# This function highlights the common promoter motifs
# Inserts a "[" before each promoter start and end, and a "]" 6 nucleotides after to create a bracketed space for each promoter motif
def AnnotateValidPromoters(DNA: str, validPromoters):
    inserts = {}  # for mapping brackets to index where it should be inserted (brackets inserted before index / see line 40, 41)
    append_end = []  # needed in case a bracket needs to be added at the end of the DNA string
    # This function inserts the brackets into the DNA sequence
    def add_insert(position, bracket):
        if position == len(DNA):
            append_end.append(bracket)
        else:
            if position not in inserts:
                inserts[position] = []
            inserts[position].append(bracket)

    for startSegment, endSegment in validPromoters:
        add_insert(startSegment, " [ ")
        add_insert(startSegment + 6, " ] ")
        add_insert(endSegment, " [ ")
        add_insert(endSegment + 6, " ] ")

    out = []
    for index, nucleotide in enumerate(DNA):
        if index in inserts:
            out.append("".join(inserts[index]))
        out.append(nucleotide)

    if append_end:
        out.append("".join(append_end))

    return "".join(out)


# These are strings that were used for testing and debugging purposes
dnaInvalid = "ACGTTTGACAXXXTATAAT" # An invalid DNA sequence because there are characters in the string that are not nucleotides
DNA = "ACGTTTGACACCGTCCCGCGCGCGCGTTCTATAAT" # A valid DNA sequence because there are 19 nucleotides between promoter motifs
DNA2 = "CAGTCAGTTTGACACGATCGGCTAGCATGTTTATAATCGATCGGGTTGACATTGCGAGCTTGACATTTTGGGGCCCGGAAAAAATTTGTATAATTCGATACGCAGTTTGACAGTTGGCAGCTAGCTTGCTATATAAT" # A valid DNA sequence with 2 promoters

if not InputValidation(DNA2):
    print("Invalid characters found. Only A,C,G,T are allowed.")
startSegments, endSegments = FindIndexOfCSeqs(DNA2)
valid = ValidateCSeqs(startSegments, endSegments)
annotated = AnnotateValidPromoters(DNA2, valid)

# This section prints out the DNA sequence with promoters highlighted
print("\nStarts:", startSegments)
print("Ends:", endSegments)
print("Valid Promoters (start_idx, end_idx):", valid)
print("Full DNA With Valid Promoters Annotated:", annotated, "\n")