from pathlib import Path
from subprocess import call
import sys
def Old_approach():
    Primary_samples = {}
    Secondary_samples = {}

    File1 = "/shares/CIBIO-Storage/CM/scratch/users/sergio.andreusanchez/Public_sample_selection/PublicDepth.tsv"; File2 = "/shares/CIBIO-Storage/CM/scratch/users/sergio.andreusanchez/DepthCoverage_Dutch.tsv"

    #Sample SGB Coverage
    for File in [File1, File2]:
     with open(File) as I:
         for line in I:
                l = line.rstrip().split()
                SGB = l[1].split("|")[-1]
                DC = float(l[2])
                if DC >= 2:
                    if SGB not in Primary_samples: Primary_samples[SGB] = []
                    Primary_samples[SGB].append(l[0])
                elif DC >0:
                    if SGB not in Secondary_samples: Secondary_samples[SGB] = []
                    Se:ndary_samples[SGB].append(l[0])

    D_prim = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/Primary/" 
    D_sec =  "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/Secondary/"

    for SGB in Primary_samples:
     Dir = D_prim + SGB + ".txt"
     with open(Dir, "a") as O:
            O.write("\n".join(Primary_samples[SGB]))
    for SGB in Secondary_samples:
        Dir = D_sec + SGB + ".txt"
        with open(Dir, "a") as O:
            O.write("\n".join(Secondary_samples[SGB]))

def New_approach(SGB_list =[], Filter=2, File1="/shares/CIBIO-Storage/CM/scratch/users/sergio.andreusanchez/Public_sample_selection/PublicDepth.tsv", File2= "/shares/CIBIO-Storage/CM/scratch/users/sergio.andreusanchez/DepthCoverage_Dutch.tsv", Location_pkl = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/All_markers/All_pkl"):
        #Dictionaries of primary and secondary samples. Per SGB, there will be a list with Paths to PKL files to use, if the depth of coverage > Filter
        Primary_samples = {};Secondary_samples = {}
        #I have symlinks of all the pkl files in Location_pkl. We will get the real path each symlink makes reference to in order to symlink in the new folder of primary and secondary samples
        PKL = {}
        print("Iterating and saving PKL locations")
        for i in Path(Location_pkl).glob("*"):
            Sample = i.stem
            Location = i.resolve()
            PKL[Sample] = Location
        #Save all SGBs in a list. This list will be used to create folders    
        SGBs = []
        #Using both files with public and private metagenomes and their read depth per bacteria (long format)
        print("Iterating through depth of covarage files. Per line, see if depth of coverage is > than filter applied")
        for File in [File1, File2]:
            with open(File) as I:
                for line in I:
                    l = line.rstrip().split()
                    SGB = l[1].split("|")[-1] #Name of the SGB. Save only the last field t_... instead of the whole name
                    if SGB not in SGB_list: continue
                    SGBs.append(SGB) #Save the SGB name
                    Sample = l[0] #Sample name is the first field
                    PKL_location = PKL[Sample] #Get the PKL location using the sample name (if it fails here, possibly trim the sample name to be whatever it is before a '_')
                    Depth = float(l[2])
                    #Now, using a specific filter (threshold of depth of coverage)
                    if Depth >= Filter: #Save primary samples above or equal than filter
                        if SGB not in Primary_samples: Primary_samples[SGB] = []
                        Primary_samples[SGB].append(PKL_location)
                    #ELIF --> Primary samples are not included in secondary samples
                    elif Depth >0: #Save secondary samples where the bacteria is present at all
                        if SGB not in Secondary_samples: Secondary_samples[SGB] = []
                        Secondary_samples[SGB].append(PKL_location)
        #Dictionary of Primary / Secondary, each key is an SGB, each entry is a location of pkl file to do symlink to
        SGB_unique = list(set(SGBs))
        if len(SGB_list)>0: SGB_unique = SGB_list
        print("Creating directories")
        for S in SGB_unique:
            #Location of SGB folder
            Location1 = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/SGB/{S}/".format(S = S)
            #location of SGB folder with specific filtering to pick primary and secondary
            Location2 = "/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/SGB/{S}/{F}".format(S = S, F= str(Filter)+"_x")
            command = "mkdir " + Location1
            command2 = "mkdir " + Location2
            if not Path(Location1).exists():
                print("Creating directory for {SGB} in : {P}".format(SGB=S, P=Location1))
                call(command, shell=True)
            if not Path(Location2).exists():
                print("Creating directory for filter in {SGB} in : {P}".format(SGB=S, P=Location2))
                call(command2, shell=True)
                

        #For primary samples. Go per SGB, get the location where PKL file should be saved. Symlink pickle file.
        for SGB in Primary_samples:
            print("Creating symlinkg for primary samples of " + SGB)
            Location =  Path("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/SGB/{S}/{R}/Primary".format(S=SGB, R=str(Filter)+"_x" ))
            if not Location.exists(): call("mkdir "+str(Location), shell=True)
            Samples = Primary_samples[SGB] #Primary samples. This is a list of Paths
            for PKL_to_link in Samples:
                PKL_to_link = Path(PKL_to_link) #In case PKL_to_link is not a Path already
                Destination = Location / PKL_to_link.name #The base name of PKL attached to the new location
                if Destination.exists(): continue #print("Symlink exists, skipping: {D}".format(D=str(Destination))) ; continue
                #print("Creating symlink {S} to {D}".format(S=str(PKL_to_link), D=str(Destination)))
                Destination.symlink_to(PKL_to_link)

        #For secondary samples. Go per SGB, get the location where PKL file should be saved. Symlink pickle file.
        for SGB in Secondary_samples:
            print("Creating symlinkg for secondary samples of " + SGB)
            Location =  Path("/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/Samples_tree/SGB/{S}/{R}/Secondary".format(S=SGB, R=str(Filter)+"_x"))
            if not Location.exists(): call("mkdir "+str(Location), shell=True)
            Samples = Secondary_samples[SGB] #Secondary samples. This is a list of Paths
            for PKL_to_link in Samples:
                PKL_to_link = Path(PKL_to_link) #In case PKL_to_link is not a Path already
                Destination = Location / PKL_to_link.name #The base name of PKL attached to the new location
                if Destination.exists(): continue #print("Symlink exists, skipping: {D}".format(D=str(Destination))) ; continue
                #print("Creating symlink {S} to {D}".format(S=str(PKL_to_link), D=str(Destination)))
                Destination.symlink_to(PKL_to_link)    

try: 
    File_sgb = sys.argv[1]
    SGBs = []
    with open(File_sgb) as O:
        for line in O:
            SGBs.append(line.rstrip())
except: SGBs = []
New_approach(SGB_list = SGBs)            
