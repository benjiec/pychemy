from pychemy.amino_acids import AMINO_ACIDS, LOSSES
from pychemy.molmass import Formula
import re


# Retreives first amino acid from peptide sequence based on those available in AMINO_ACIDS
def get_next_aa(sequence):
  keep = []
  for key in AMINO_ACIDS:
    if sequence.startswith(key.mod):
      keep.append(key)
  
  if len(keep) == 1:
    return keep[0]
  elif len(keep) == 0:
    raise Exception('No matching residue in amino acid set')
  else:
    mod_set = [k for k in keep if k.mod != k.residue]
    if len(mod_set) == 1:
      return mod_set[0]
    else:
      raise Exception('Conflicting modifications in amino acid set')

# Parse string containing sequence into residues
def residues_from_sequence(sequence):
  count = 0
  max_iter = len(sequence)
  residues = []
  while len(sequence) > 0 and count < max_iter + 2:
    count += 1
    aa = get_next_aa(sequence)
    residues.append(aa.mod)
    sequence = sequence.replace(aa.mod, '', 1)
  return residues

# Check that a provided sequence is valid
def verify_sequence(sequence):
  assert (type(sequence) is not str), 'Sequence must be an array of amino acids, which can be produced from residues_from_sequence'

  for r in sequence:
    assert (r in [key.mod for key in AMINO_ACIDS]), 'Sequence must contain valid amino acids as defined in AMINO_ACIDS'

  return True

# Calculate mass of peptide based on chemical formula entered as 'ACDE'
def mass_from_sequence(sequence='', nterm='N-term', cterm='C-term'):
  seq_residues = residues_from_sequence(sequence)
  mass = 0

  # Include N-terminal group mass
  mass += AMINO_ACIDS[nterm].mass

  # Include each residue mass
  for r in seq_residues:
    mass += AMINO_ACIDS[r].mass
      
  # Include C-terminal group mass
  mass += AMINO_ACIDS[cterm].mass

  return mass

def get_precursors(seq, 
                   min_mz=400, 
                   max_mz=2000, 
                   min_charge=1, 
                   max_charge=5):
  m = mass_from_sequence(seq)
  prec = [{
           'charge_state': c, 
           'mz': (m + c * Formula('H').mass) / c
          } for c in xrange(min_charge, max_charge + 1)]
  return filter(lambda x: min_mz <= x['mz'] <= max_mz, prec)

def update_loss(loss, new_loss):
  for idx,l in enumerate(loss):
    if l[1].formula == LOSSES[new_loss].formula:
      loss[idx][0] += 1
      return loss
  return loss + [[1, LOSSES[new_loss]]]


def loss_name(loss=[]):
  out = ''
  for l in loss:
    if l[0] == 1:
      out += ' - ' + l[1].display_name
    else:
      out += ' - ' + str(l[0]) + ' ' + l[1].display_name
  return out


def predict_positive_fragments_from_sequence(sequence=[], nterm='N-term', cterm='C-term',
                                             max_charge_state=1):
  out = []

  fragments = predict_fragments_from_sequence(sequence=sequence, nterm=nterm, cterm=cterm)

  for item in fragments:
    out.append(( round(item[0] + Formula('H').isotope.mass, 5) , item[1] + ' (+1)'))
    
    for cs in range(2,max_charge_state + 1):
      if item[0] > 200 * cs:
        out.append(( round((item[0] + cs * Formula('H').isotope.mass) / cs, 5), item[1] + ' (+' +  str(cs) + ')'))
  
  out = sorted(list(set(out)), key = lambda x:x[0])
  return out


def predict_negative_fragments_from_sequence(sequence=[], nterm='N-term', cterm='C-term',
                                             max_charge_state=1):
  out = []

  fragments = predict_fragments_from_sequence(sequence=sequence, nterm=nterm, cterm=cterm)

  for item in fragments:
    out.append(( round(item[0] - Formula('H').isotope.mass, 5), item[1] + ' (-1)'))
    
    for cs in range(2,max_charge_state + 1):
      if item[0] > 200 * cs:
        out.append(( round((item[0] - cs * Formula('H').isotope.mass) / cs, 5), item[1] + ' (-' +  str(cs) + ')'))
  
  out = sorted(list(set(out)), key=lambda x:x[0])
  return out


# Predict fragments from peptide sequence    
def predict_fragments_from_sequence(sequence=[], nterm='N-term', cterm='C-term'):
  
  verify_sequence(sequence)

  out = []

  # b-ions
  nterm_mass = AMINO_ACIDS[nterm].mass
  out += append_residue(mass=nterm_mass, 
                        pos=1, 
                        ion_type='b_', 
                        seq=sequence, 
                        term=cterm, 
                        loss=[])
  for l in AMINO_ACIDS[nterm].losses:
    out += append_residue(mass=nterm_mass - l.mass, 
                          pos=1, 
                          ion_type='b_', 
                          seq=sequence, 
                          term=cterm, 
                          loss=update_loss([],l.formula))

  # a-ions
  nterm_mass = AMINO_ACIDS[nterm].mass - Formula('CO').isotope.mass 
  out += append_residue(mass=nterm_mass, 
                        pos=1, 
                        ion_type='a_', 
                        seq=sequence, 
                        term=cterm, 
                        loss=[])
  for l in AMINO_ACIDS[nterm].losses:
    out += append_residue(mass=nterm_mass - l.mass, 
                          pos=1, 
                          ion_type='a_', 
                          seq=sequence, 
                          term=cterm, 
                          loss=update_loss([],l.formula))

  # y-ions
  cterm_mass = AMINO_ACIDS[cterm].mass
  out += append_residue(mass=cterm_mass, 
                        pos=1, 
                        ion_type='y_', 
                        seq=sequence[::-1], 
                        term=nterm, 
                        loss=[])
  for l in AMINO_ACIDS[cterm].losses:
    out += append_residue(mass=cterm_mass - l.mass, 
                          pos=1, 
                          ion_type='y_', 
                          seq=sequence[::-1], 
                          term=nterm, 
                          loss=update_loss([],l.formula))
      
  out = sorted(list(set(out)), key=lambda x:x[0])
  return out


def append_residue(mass, pos, ion_type='', seq=[], term='', loss=[]):
  out = []

  if pos == 1 and not seq:
    raise Exception("append_residue encountered empty sequence")

  if len(seq) == 0:
    # No residues left to place, place terminus
    r = AMINO_ACIDS[term]

    if ion_type ==  'a_':
      # Special case: M - CO
      loss = update_loss(loss, 'CO')  
    
    out += [(round(mass + r.mass, 5), 'M' + loss_name(loss))]
    for l in r.losses: 
      out += [(round(mass + r.mass - l.mass, 5), 'M' + loss_name(update_loss(loss, l.formula)))]

  else:
    # Place each loss form of next residue
    r = AMINO_ACIDS[seq[0]]

    out += [(round(mass + r.mass, 5), ion_type + str(pos) + loss_name(loss))]          
    out += append_residue(mass=mass + r.mass, 
                          pos=pos + 1, 
                          ion_type=ion_type, 
                          seq=seq[1:], 
                          term=term, 
                          loss=loss)
    for l in r.losses:
      curr_loss = update_loss(loss, l.formula)
      out += [(round(mass + r.mass - l.mass, 5), ion_type + str(pos) + loss_name(curr_loss))]
      out += append_residue(mass=mass + r.mass - l.mass, 
                            pos=pos + 1, 
                            ion_type=ion_type, 
                            seq=seq[1:], 
                            term=term, 
                            loss=curr_loss)

  return out
