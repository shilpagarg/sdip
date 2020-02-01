import sys
import pysam
import pandas as pd
import argparse
import networkx as nx
from collections import Counter


def form_graph(df, fasta_file, identity_cutoff, max_similarity):
    #Filter out duplicate contigs
    identical_contigs = df.loc[(df["perID_by_all"] >= max_similarity) & (df["query_span"] / df["query_length"] > 0.95) & (df["reference_span"] / df["reference_length"] > 0.95) & (df["query_name"] < df["reference_name"]), ["query_name", "reference_name"]]
    duplicates = set([row["query_name"] for index, row in identical_contigs.iterrows()])
    contig_names = [contig for contig in fasta_file.references if contig not in duplicates]
    
    similar_contigs = df.loc[(df["perID_by_all"] >= identity_cutoff) & (df["query_span"] / df["query_length"] > 0.5) & (df["reference_span"] / df["reference_length"] > 0.5) & ~(df["query_name"].isin(duplicates)) & ~(df["reference_name"].isin(duplicates)), ["query_name", "reference_name"]]
    
    #Build graph
    G = nx.Graph()
    for contig in contig_names:
        G.add_node(contig, length = fasta_file.get_reference_length(contig))
    for index, row in similar_contigs.iterrows():
        G.add_edge(row["query_name"], row["reference_name"])
    return G, duplicates

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('contigs', type=str, help = 'Contigs in FASTA format')
    parser.add_argument('tbl', type=str, help = 'Table from samIdentities.py')
    parser.add_argument('region', type=int, help = 'Region number for FASTA sequence headers')
    parser.add_argument('--min_similarity', '-m', type=float, default = 98.0, help = 'Minimum similarity cutoff for contig pairing (default: 98.0)')
    parser.add_argument('--max_similarity', '-x', type=float, default = 99.9, help = 'Contigs more similar than this value are merged (default: 99.9)')
    parser.add_argument('--haploid', action='store_true', help = 'Organism is haploid')
    args = parser.parse_args()

    #Read contig names
    fasta_file = pysam.FastaFile(args.contigs)
    
    #Read table
    df = pd.read_csv(args.tbl, sep="\t")
    df["query_span"] = df["query_end"] - df["query_start"]
    df["reference_span"] = df["reference_end"] - df["reference_start"]

    if args.haploid:
        identical_contigs = df.loc[(df["perID_by_all"] >= args.max_similarity) & (df["query_span"] / df["query_length"] > 0.95) & (df["reference_span"] / df["reference_length"] > 0.95) & (df["query_name"] < df["reference_name"]), ["query_name", "reference_name"]]
        duplicates = set([row["query_name"] for index, row in identical_contigs.iterrows()])
        contig_names = [contig for contig in fasta_file.references if contig not in duplicates]
        for d in duplicates:
            print("DUPLICATES\t%s" % (d), file=sys.stderr)
        p_index = 0
        for contig in contig_names:
            p_index += 1
            #Put single contigs into haplotype0
            print(">r%s_paralog%s" % (args.region, p_index))
            print(fasta_file.fetch(reference = contig))
    else:
        start = args.min_similarity
        end = 100.0
        step = 0.05
        steps = int((end - start) / step)
        max_pairs = -1
        best_cutoff = -1
        best_components = None
        best_duplicates = None
        for cutoff in [start + i * step for i in range(steps)]:
            graph, duplicates = form_graph(df, fasta_file, cutoff, args.max_similarity)
            num_pairs = len([len(component) for component in nx.connected_components(graph) if len(component) == 2])
            if num_pairs > max_pairs:
                max_pairs = num_pairs
                best_cutoff = cutoff
                best_components = list(nx.connected_components(graph))
                best_duplicates = duplicates
            print("CUTOFF\t%s\t%s" % (cutoff, [(len(component), graph.nodes[list(component)[0]]["length"]) for component in nx.connected_components(graph)]), file=sys.stderr)
        
        c = Counter()
        for l in [len(component) for component in best_components]:
            c[l] += 1
        for l, count in c.most_common():
            print("COMPONENTS\t%d\t%d" % (l, count), file=sys.stderr)
        for d in best_duplicates:
            print("DUPLICATES\t%s" % (d), file=sys.stderr)
        
        p_index = 0
        for component in best_components:
            p_index += 1
            if len(component) == 1:
                #Put single contigs into haplotype0
                contig = list(component)[0]
                print(">r%s_paralog%s_haplotype0" % (args.region, p_index))
                print(fasta_file.fetch(reference = contig))
            elif len(component) == 2:
                #Put paired contigs into haplotype1 and haplotype2
                for c_index, contig in enumerate(component):
                    print(">r%s_paralog%s_haplotype%s" % (args.region, p_index, c_index+1))
                    print(fasta_file.fetch(reference = contig))
            else:
                for c_index, contig in enumerate(component):
                    #Every second contig goes into a separate paralog
                    p_index2 = c_index // 2
                    #Put last single contig into haplotype0
                    if (c_index + 1) == len(component) and len(component) % 2 == 1:
                        print(">r%s_paralog%s_haplotype0" % (args.region, p_index + p_index2))
                        print(fasta_file.fetch(reference = contig))
                    else:
                        print(">r%s_paralog%s_haplotype%s" % (args.region, p_index + p_index2, (c_index % 2) + 1))
                        print(fasta_file.fetch(reference = contig))
                p_index+=p_index2

if __name__ == '__main__':
    main()
