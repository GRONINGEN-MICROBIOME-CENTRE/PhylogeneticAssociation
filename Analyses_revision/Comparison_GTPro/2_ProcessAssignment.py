from pathlib import Path

Info_ass = "Assignment/Assignment.tsv"
Info = dict()
with open(Info_ass, "r") as File:
    for line in File:
        if line[0] == "#": continue
        l = line.rstrip().split()
        Genome = l[0].split(".")[0]
        SGB = l[1].split("|")[-1]
        Dist = round(float(SGB.split(":")[-1]),2)
        SGB = SGB.split(":")[0]
        FASTA = "../../" + l[0] + ".fasta"
        Distance_file = "../../../../" + Genome + ".distance.RData"
        if not Path(FASTA).exists():
            exit('Error: ' + FASTA + ' does not exist')
        Info[Distance_file] = { 'Bug' : Genome, 'SGB' : SGB, 'Distance_ref' : str(Dist)   } 

Symlink_loc = Path("/scratch/p300317/Strains_sergio/Analysis/1_AssignTaxonomy/Assignment/Genomes_assigned")
for Location in Info:
    I = Info[Location]
    Location = Path(Location)
    NewLocation = Symlink_loc / (I["SGB"] + ".RData")
    Distance = I["Distance_ref"]
    if float(Dist) > 0.5: exit(I)
    if not NewLocation.exists():
        NewLocation.symlink_to(Location)
