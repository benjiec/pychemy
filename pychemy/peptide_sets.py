# Tools for grouping proteins into non-conflicting sets for
# multiplexed screening of homologues

import re

AA = re.compile('[ACDEFGHILMNPQSTVWY]+[RK]')
AA_check = re.compile('[ACDEFGHILKMNPQRSTVWY]')

# Check that all entries in sequence represent valid amino acids
def check_AA(test_string):
  if len(AA_check.findall(test_string)) == len(test_string):
    return True
  else:
    raise Exception("Invalid amino acid in seq: " + test_string)

# Check that all proteins in list contain only valid amino acids
# input: 
#        proteins: list of protein sequences as strings ['protein_seq', ...]
def check_proteins(proteins):
  for p in proteins:
    check_AA(p)
  return True

# Creates list of peptides from protein sequence
# input:
#        protein_seq: string containing protein sequence of amino acids
# output:
#        list of strings containing peptide sequences from protein sequence ['peptide_seq', ...]
def get_peptides(protein_seq):
  check_AA(protein_seq)
  return AA.findall(protein_seq)

# Creates list of all peptides from protein set
# inputs:
#        proteins: list of protein sequences as strings ['protein_seq', ...]
# output:
#        list of strings containing peptide sequences from all protein sequences ['peptide_seq', ...] 
def all_peptides(proteins):
  return [p for peptides in [get_peptides(ps) for ps in proteins] for p in peptides]

# Creates list of all distinct peptides from protein set
# inputs:
#        proteins: list of protein sequences as strings ['protein_seq', ...]
# output:
#        list of strings containing distinct peptide sequences from all protein sequences ['peptide_seq', ...] 
def distinct_peptides(proteins):
  return list(set(all_peptides(proteins)))


# Creates list of unique identifiers for each protein from set
# inputs:
#        proteins: list of protein sequences as strings ['protein_seq', ...]
# output:
#        list of lists containing protein sequence and a unique identifying peptide sequence [['protein_seq', 'peptide_seq'], ...]       
def unique_identifiers(proteins):
    all_pep = [[ps[0], pep] for ps in enumerate(proteins) for pep in get_peptides(ps[1])]
    dp = distinct_peptides(proteins)
    ui = []
    for idx,seq in enumerate(dp):
      idx = [item[0] for item in all_pep if item[1] == seq]
      if len(idx) == 1:
        ui.append([proteins[idx[0]], seq])
    return ui

# Creates list of unique identifier peptides grouped by protein
# input: 
#        ui: a set of unique identifiers such as that produced by unique_identifiers([proteins]) [['protein_seq', 'peptide_seq'],.....]
# output: 
#        list of unique identifier peptides grouped by their protein [['protein_seq', ['peptide_seq', 'peptide_seq', ...]],.....]
def peptides_per_protein(ui):
  proteins = list(set([item[0] for item in ui]))
  return [[p, [item[1] for item in ui if item[0] == p]] for p in proteins]

# Counts the minimum number of unique identifier peptides across all proteins in a set
# input:
#        proteins: list of protein sequences as strings ['protein_seq', ...]
# output:
#        minimum number of unique identifier peptides across all proteins in a set
def min_ui_count(proteins):
  temp = []
  for p in peptides_per_protein(unique_identifiers(proteins)):
    temp.append(len(p[1]))
  if len(proteins) > len(temp):
    return 0
  else:
    return min(temp)

# Creates balanced sets of proteins based on unique peptides subject to contraints
# inputs:
#        proteins:     list of protein sequences as strings ['protein_seq', ...]
#        max_set_size: maximum number of proteins permitted in a set
#        unique:       minimum number of unique peptides required per protein
# output:
#        sets of protein sequences that meet the set size and uniqueness requirements [['protein_seq', ...], ...]
def balanced_sets(proteins, max_set_size = 10, unique = 1):
  check_proteins(proteins)
  sets = []
  for p in proteins:
    placed = False
    for idx,s in enumerate(sets):
      if len(s) < max_set_size and not placed:
        if min_ui_count(s + [p]) >= unique:
          placed = True
          sets[idx] += [p]   
    if not placed:
     sets.append([p])
  return sets

