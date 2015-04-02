from elements import ELEMENTS
import networkx as nx
from collections import OrderedDict
import openbabel 
import matplotlib.pyplot as plt

class chem_graph():

  def __init__(self, OBMol = None, graph = None, is_fragment = False, parent = None):
    G = None
    if OBMol:
      G = chem_graph.graph_from_OBMol(OBMol)      
    elif graph:
      G = nx.Graph(graph)
    
    self.G = nx.Graph(G)

    if is_fragment:
      self.parent = parent
    else:
      self.parent = self

  def gen_frag(self, num = 1):
    if num == 0:
      return [self]
    else:
      frag = [self]
      for item in self.G.edges():
        G = nx.Graph(self.G)
        G.remove_edge(item[0],item[1])
        while G.number_of_nodes() > 0:
          include = [i for i in nx.dfs_preorder_nodes(G,G.nodes()[0])]
          exclude = [e for e in G.nodes() if e not in include]
          G2 = nx.Graph(G)
          for n in exclude:
            G2.remove_node(n)
          frag += chem_graph(graph=G2, is_fragment = True, parent = self.parent).gen_frag(num-1)
          for n in include:
            G.remove_node(n)
          frag += chem_graph(graph=G, is_fragment = True, parent = self.parent).gen_frag(num-1)
      return frag

  def chem_formula(self):
    atoms = dict()
    for n in self.G.nodes():
      a = self.G.node[n]['atom'].symbol
      if a in atoms:
        atoms[a] += 1
      else:
        atoms[a] = 1

    atoms = OrderedDict(sorted(atoms.items(), key = lambda t:t[0]))

    cf = ''
    for key in atoms:
       cf += key + ' ' + str(atoms[key]) + ' '
    self.atoms = atoms
    return cf[:-1]

  def mass(self):
    cf = self.chem_formula()
    m = 0
    for key in self.atoms:
      m += self.atoms[key] * ELEMENTS[key].exactmass
    return m

  @staticmethod
  def graph_from_OBMol(mol):
    G = nx.Graph()

    # Generate 2d positions
    openbabel.OBOp.FindType("gen2d").Do(mol)
    
    for atom in openbabel.OBMolAtomIter(mol):
      idx = atom.GetIdx()

      G.add_node(idx)
      G.node[idx]['atom'] = ELEMENTS[atom.GetAtomicNum()]

      # Save 2d positions for later plotting
      G.node[idx]['X'] = atom.GetX()
      G.node[idx]['Y'] = atom.GetY()

    for bond in openbabel.OBMolBondIter(mol):
      atom1 = bond.GetBeginAtom().GetIdx()
      atom2 = bond.GetEndAtom().GetIdx()
      G.add_edge(atom1, atom2)
      G[atom1][atom2]['order'] = bond.GetBO()  
  
    return G
  
  def draw(self, node_color = (1,1,1), edge_color = [0,0,0], node_size = 700, font_color = [0,0,0]):
    node_colors = [node_color for i in range(self.G.number_of_nodes())]
    edge_colors = [edge_color for i in range(self.G.number_of_edges())]
    widths = [2 for i in range(self.G.number_of_edges())]
    pos_dict = dict()
    label_dict = dict()
    for n in self.G.nodes():
      pos_dict[n] = [self.G.node[n]['X'], self.G.node[n]['Y']]
      label_dict[n] = self.G.node[n]['atom'].symbol
    
    # Draw single bonds
    nx.draw(self.G,
            pos = pos_dict,
            labels = label_dict,
            with_labels = True,
            font_size = 24,
            font_color = font_color,
            node_size = node_size,
            node_color = node_colors,
            edge_color = edge_colors,
            linewidths = 0,
            width = widths,
            )

    # Draw double bonds
    edge_list = [edge for edge in self.G.edges() if self.G[edge[0]][edge[1]]['order'] == 2]
    nx.draw(self.G,
            pos = pos_dict,
            labels = label_dict,
            with_labels = True,
            font_size = 24,
            font_color = font_color,
            node_size = node_size,
            node_color = node_colors,
            edge_color = [edge_color for e in edge_list],
            linewidths = 0,
            width = [6 for i in widths],
            edgelist = edge_list
            )

    nx.draw(self.G,
            pos = pos_dict,
            labels = label_dict,
            with_labels = True,
            font_size = 24,
            font_color = font_color,
            node_size = node_size,
            node_color = node_colors,
            edge_color = [[1,1,1] for e in edge_list],
            linewidths = 0,
            width = [2 for i in widths],
            edgelist = edge_list
            )
    
    # Draw triple bonds
    edge_list = [edge for edge in self.G.edges() if self.G[edge[0]][edge[1]]['order'] == 3]
    nx.draw(self.G,
            pos = pos_dict,
            labels = label_dict,
            with_labels = True,
            font_size = 24,
            font_color = font_color,
            node_size = node_size,
            node_color = node_colors,
            edge_color = [edge_color for e in edge_list],
            linewidths = 0,
            width = [10 for i in widths],
            edgelist = edge_list
            )

    nx.draw(self.G,
            pos = pos_dict,
            labels = label_dict,
            with_labels = True,
            font_size = 24,
            font_color = font_color,
            node_size = node_size,
            node_color = node_colors,
            edge_color = [[1,1,1] for e in edge_list],
            linewidths = 0,
            width = [6 for i in widths],
            edgelist = edge_list
            )
     
    nx.draw(self.G,
            pos = pos_dict,
            labels = label_dict,
            with_labels = True,
            font_size = 24,
            font_color = font_color,
            node_size = node_size,
            node_color = node_colors,
            edge_color = [edge_color for e in edge_list],
            linewidths = 0,
            width = [2 for i in widths],
            edgelist = edge_list
            )

  def png(self, fn):
    if self.parent:
      w = 0.95
      self.parent.draw(edge_color = [w,w,w],
                font_color = [w,w,w], 
                node_size = 600,
                )
    self.draw()
    plt.savefig(fn+'.png')
    plt.close()


class chem_structure():
  def __init__(self, inchi = None):
    self.inchi = inchi
    self.mol = openbabel.OBMol()
    if inchi:
      obConversion = openbabel.OBConversion()
      obConversion.SetInAndOutFormats("inchi", "mdl")
      obConversion.ReadString(self.mol, inchi)
      self.mol.AddHydrogens()
 
  def chem_formula(self):
    return self.mol.GetSpacedFormula()

  def mass(self):
    return self.mol.GetExactMass()

  def fragments(self, bonds = 1):
    CG = chem_graph(OBMol = self.mol)
    u = dict()
    for item in CG.gen_frag(num = bonds):
      nodes = tuple(sorted([n for n in item.G.nodes()]))
      if nodes not in u.keys():
        u[nodes] = item
    return [u[key] for key in u]  
