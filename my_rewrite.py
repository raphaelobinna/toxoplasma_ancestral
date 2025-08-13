import pandas as pd
from Bio import SeqIO, Phylo

def read_lines(file_path):
    lines =[]
    with open(file_path, 'r') as f:
        for line in f:
            line =line.strip()
            if line and not line.startswith('#'):
                lines.append(line)

def read_csv_file(file_path):
    ancestral_data = pd.read_csv(file_path, sep='\t', comment='#')
    print(f'Columns are {list(ancestral_data.columns)}')
    print(f"Unique Nodes length {ancestral_data['Node'].nunique()}")
    return ancestral_data


def read_sequences(file_path):
    sequences={}
    headers={}

    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
        headers[record.id] = record.description
    print(f'The headers of sequences {headers}')
    return sequences, headers

def read_trees(file_path):
    tree = Phylo.read(file_path, "newick")
    print(f'treee termiansls = {tree.get_terminals()}')

    return tree

def get_tree_structure(tree):

    edges =[]

    def looper(clade, parent_name=None):
        if clade.name is None:
            clade.name = f'New_{id(clade)}'
        if parent_name:
            edges.append((parent_name, clade.name))

        for child in clade.clades:
            looper(child, clade.name)
        

    looper(tree.root)
    return edges

def compare_parent_child(ancestor, branches):
    mutations = []
    states = {(row.Node, row.Site):row for _, row in ancestor.iterrows()}
    unique_sites = ancestor.Site.unique()

    for parent, child in branches:
        for site in unique_sites:
            if (parent, site) in states and (child, site) in states:
                p_state = states[(parent, site)]
                c_state = states[(child, site)]

                if p_state.State != c_state.State:
                    mutations.append({
                        "Branch": f"{parent}->{child}",
                        "Site": f"{site}",
                        "Change": f"{p_state.State}{site}{c_state.State}",
                        "Parent_p_C": p_state.p_C,
                        "Parent_p_G": p_state.p_G,
                        "Child_p_C": c_state.p_C,
                        "Child_p_G": c_state.p_G
                    })

    return pd.DataFrame(mutations)

            

def main():
    """Main analysis pipeline"""
    # File paths - adjust these to your files
    ancestral_file = 'dhfr_codon_aln.fasta.state'
    tree_file = 'dhfr_codon_aln.fasta.treefile'
    sequence_file = 'dhfr_codon_aln.fasta'


    try:

        ancestral_data_read = read_csv_file(ancestral_file)
        tree_data = read_trees(tree_file)
        sequence_data = read_sequences(sequence_file)

        all_branches = get_tree_structure(tree_data)
        print(f'Number of branches {len(all_branches)}')

        all_mutations = compare_parent_child(ancestral_data_read, all_branches)
        print(f"data head {all_mutations.head()}")
        all_mutations.to_csv("mutations_simple.csv", index=False)




    except FileNotFoundError as e:
        print(f"File not found: {e}")
        print("Please check your file paths and make sure all files exist")
    except Exception as e:
        print(f"Error during analysis: {e}")
        return None

if __name__ == "__main__":
    # Run the analysis
    results = main()