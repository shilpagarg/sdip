import subprocess
from collections import defaultdict

with open(snakemake.input["contigs"], 'r') as contigs_file:
    content = contigs_file.read()
if len(content) < 1:
    subprocess.run(["touch", snakemake.output["contigs"]])
else:
    #Extract contigs into separate files
    subprocess.run(["mkdir", "-p", snakemake.params["wd"]])
    num = 0
    contig_name = ""
    sequence = ""
    for line in content.split("\n"):
        if line.startswith(">"):
            if contig_name != "" and sequence != "":
                num += 1
                with open(snakemake.params["wd"] + "contig%d.fa" % (num), "w") as contig_file:
                    print(">" + contig_name, file=contig_file)
                    print(sequence, file=contig_file)
            contig_name = line.strip()[1:].split()[0]
            sequence = ""
        else:
            sequence += line.strip()
    if contig_name != "" and sequence != "":
        num += 1
        with open(snakemake.params["wd"] + "contig%d.fa" % (num), "w") as contig_file:
            print(">" + contig_name, file=contig_file)
            print(sequence, file=contig_file)
    num_haps = num
    
    #Store reads for each contig
    with open(snakemake.input["reads"], 'r') as read_file:
        reads = defaultdict(list)
        for line in read_file:
            fields = line.strip().split()
            num = int(fields[0])
            assert num <= num_haps
            reads[num].append(fields[1])
    
    #Add contained reads
    with open(snakemake.input["contained_reads"], 'r') as contained_file:
        for num in range(1, num_haps + 1):
            contained_file.seek(0)
            for line in contained_file:
                fields = line.strip().split()
                if fields[1] in reads[num]:
                    reads[num].append(fields[0])
    
    #Polish each contig
    for num in range(1, num_haps + 1):
        reads_path = snakemake.params["wd"] + "reads%d.txt" % (num)
        with open(reads_path, 'w') as read_file:
            for r in reads[num]:
                print(r, file=read_file)
        fasta_path = snakemake.params["wd"] + "reads%d.fa" % (num)
        with open(fasta_path, "w") as outfile:
            subprocess.run(["samtools", "faidx", "-r", reads_path, snakemake.input["sequences"]], stdout = outfile)
        contig_path = snakemake.params["wd"] + "contig%d.fa" % (num)
        sam_path = snakemake.params["wd"] + "aligned%d.sam" % (num)
        with open(sam_path, "w") as outfile:
            subprocess.run(["minimap2", "-ax", "map-pb", "-t4", contig_path, fasta_path], stdout = outfile)
        consensus_path = snakemake.params["wd"] + "consensus%d.fa" % (num)
        with open(consensus_path, "w") as outfile:
            subprocess.run(["racon", "--no-trimming", "-t", str(snakemake.threads), "-u", fasta_path, sam_path, contig_path], stdout = outfile)
    
    #Concat polished contigs
    with open(snakemake.output["contigs"], "w") as outfile:
        subprocess.run(["cat"] + [(snakemake.params["wd"] + "consensus%d.fa" % (num)) for num in range(1, num_haps + 1)], stdout = outfile)
