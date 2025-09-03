def InputValidation(DNA: str) -> bool:
    for c in DNA:
        if c not in ("A", "C", "G", "T"):
            return False
    return True

def FindIndexOfCSeqs(DNA: str):
    start_motif = "TTGACA"  # start of promoter (-35)
    end_motif   = "TATAAT"  # end of promoter (-10)
    indexOfPotentialStartPromoters = []
    indexOfPotentialEndPromoters = []
    dnaLength = len(DNA)
    for i in range(0, dnaLength - 5):
        segment = DNA[i:i+6]
        if segment == start_motif:
            indexOfPotentialStartPromoters.append(i)
        if segment == end_motif:
            indexOfPotentialEndPromoters.append(i)
    return indexOfPotentialStartPromoters, indexOfPotentialEndPromoters

def ValidateCSeqs(starts, ends):
    # Valid if end comes after start and there are 16–19 characters between them,
    validPromoters = []
    for s in starts:
        for e in ends:
            spacer = e - (s + 6)
            if 16 <= spacer <= 19:
                validPromoters.append((s, e))
    return validPromoters

def AnnotateValidPromoters(DNA: str, validPromoters):
    # Insert " [ " before each promoter start and end, and " ] " 6 neucleotides after
    inserts = {}  # for mapping brackets to index where it should be inserted (brackets inserted before index / see line 40, 41)
    append_end = []  # needed in case bracket needs to be added at end of DNA string
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


# Testing / Debugging
dnaInvalid = "ACGTTTGACAXXXTATAAT"
DNA = "ACGTTTGACACCGTCCCGCGCGCGCGTTCTATAAT" # 19 bps between segments
DNA2 = "CAGTCAGTTTGACACGATCGGCTAGCATGTTTATAATCGATCGGGTTGACATTGCGAGCTTGACATTTTGGGGCCCGGAAAAAATTTGTATAATTCGATACGCAGTTTGACAGTTGGCAGCTAGCTTGCTATATAAT" # Should be 2 valid motifs

if not InputValidation(DNA2):
    print("Invalid characters found. Only A,C,G,T are allowed.")
startSegments, endSegments = FindIndexOfCSeqs(DNA2)
valid = ValidateCSeqs(startSegments, endSegments)
annotated = AnnotateValidPromoters(DNA2, valid)

print("\nStarts:", startSegments)
print("Ends:", endSegments)
print("Valid pairs (start_idx, end_idx):", valid)
print("Valid consensus sequences:", annotated, "\n")