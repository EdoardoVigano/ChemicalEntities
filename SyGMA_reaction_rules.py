from rdkit import Chem

'''Note:
Every time you add something here you have to add the same changes in thesefile.py: 

 1.alerts_index
 2.calculator 
 3.alerts
'''
def REACTION_N_dealkylation_1(molecule: Chem.rdchem.Mol):
	# 'N-demethylation_(R-NHCH3)'
	substructure = '[*;!c:1][NH1;X3:2][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_2(molecule: Chem.rdchem.Mol):
	# 'N-demethylation_(c-NHCH3)'
	substructure = '[c:1][NH1;X3:2][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_3(molecule: Chem.rdchem.Mol):
	# 'N-demethylation_(R-N(CH3)2)'
	substructure = '[*;!c:1][NH0;X3:2]([CH3])[CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_4(molecule: Chem.rdchem.Mol):
	# 'N-demethylation_(c-N(CH3)2)'
	substructure = '[c:1][NH0;X3:2]([CH3])[CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_5(molecule: Chem.rdchem.Mol):
	# 'N-demethylation_(R-N(CR)CH3)'
	substructure = '[*;!$([CH3]):1][NH0;X3:2]([CH3])[#6;!$([CH3]):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_6(molecule: Chem.rdchem.Mol):
	# 'N-demethylation_(nCH3)'
	substructure = '[n:1][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_7(molecule: Chem.rdchem.Mol):
	# 'N-depropylation'
	substructure = '[N;X3:2][CH1]([CH3])[CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_8(molecule: Chem.rdchem.Mol):
	# 'secondary_N-depropylation'
	substructure = '[NH1;X3:2][CH1]([CH3])[CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_9(molecule: Chem.rdchem.Mol):
	# 'tertiary_N-depropylation'
	substructure = '[NH0;X3:2][CH1]([CH3])[CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_10(molecule: Chem.rdchem.Mol):
	# 'N-deglycosidation'
	substructure = '[NX3:2][C:3]1[O:4][C:5][C:6][C:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_11(molecule: Chem.rdchem.Mol):
	# 'n-deglycosidation'
	substructure = '[n:2][C:3]1[O:4][C:5][C:6][C:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_12(molecule: Chem.rdchem.Mol):
	# 'N-deformylation'
	substructure = '[NX3:2][CX3;H1]=O'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_13(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(piperazine)'
	substructure = '[*;!C,!X4:1][N;X3:2]1[C:3][C:4][N;X3:5][CH2][CH2]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_14(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(morpholine)'
	substructure = '[N;X3:2]1[C:3][C:4][O:5][CH2][CH2]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_15(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(R-NHCH2-alkyl)'
	substructure = '[*;!c:1][NH1;X3:2]!@[CH2:3][#6:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_16(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(c-NHCH2-alkyl)'
	substructure = '[c:1][NH1;X3:2]!@[CH2:3][#6:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_17(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(tertiaryN-CH2-alkyl)'
	substructure = '[NH0;X3:2]!@[C;X4;H2:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_18(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(quarternary_N)'
	substructure = '[#6:1][N+;X4:2]([#6:3])([CH3:4])!@[#6;H1,H2:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_19(molecule: Chem.rdchem.Mol):
	# 'N-dealkylation_(nCH2)'
	substructure = '[n:1][CH2:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_20(molecule: Chem.rdchem.Mol):
	# 'tertiary_N-dealkylation2'
	substructure = '[NH0;X3:2]!@[C;X4;H1:4][c:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_dealkylation_21(molecule: Chem.rdchem.Mol):
	# 'tertiary_N-dealkylation'
	substructure = '[#6:1][N:2]([#6:3])[c:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_1(molecule: Chem.rdchem.Mol):
	# 'het-O-demethylation'
	substructure = '[#6:1]!@[N;R;X3:2]([CH2:3])[CH2:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_2(molecule: Chem.rdchem.Mol):
	# 'O-dealkylation_(aliphatic)'
	substructure = '[#6;!$(C=O):1][O:2][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_3(molecule: Chem.rdchem.Mol):
	# 'O-dealkylation_(aromatic)'
	substructure = '[*;!#6;!$(*=O):1][O:2][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_4(molecule: Chem.rdchem.Mol):
	# 'O-deglycosidation'
	substructure = '[C;!$(C(O)~[!#6]);!$([CH3]):1][O;!$(O1CC1):2][C;X4;!$(C(O)~[!#6]);H1,H2:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_5(molecule: Chem.rdchem.Mol):
	# 'O-dealkylation_(methylenedioxyphenyl)a'
	substructure = '[c:1][O:2][C;X4;!$(C(O)~[!#6]);H1,H2:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_6(molecule: Chem.rdchem.Mol):
	# 'O-dealkylation_(methylenedioxyphenyl)b'
	substructure = '[#6;!$([CH3]);!$(C=O):1][O:2][C:3]1[O:4][C:5][C:6][C:7][C:8]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_7(molecule: Chem.rdchem.Mol):
	# '-'
	substructure = '[O:1]1[c:2]2[c:3][c:4][c:5][c:6][c:7]2[O:8][CH2]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_demethylation_8(molecule: Chem.rdchem.Mol):
	# '-'
	substructure = '[O:1]1[c:2]2[c:3][c:4][c:5][c:6][c:7]2[O:8][CH2:9]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_dealkylation_1(molecule: Chem.rdchem.Mol):
	# 'S-dealkylation_c-SCH2-R'
	substructure = '[c:1][S:2][CH2:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_1(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(general)'
	substructure = '[cH1:1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_2(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(para_to_carbon)'
	substructure = '[#6:1]~[a:2]1[a:3][a:4][cH1:5][a:6][a:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_3(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(para_to_nitrogen)'
	substructure = '[#7:1]~[a:2]1[a:3][a:4][cH1:5][a:6][a:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_4(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(para_to_oxygen)'
	substructure = '[#8:1]~[a:2]1[a:3][a:4][cH1:5][a:6][a:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_5(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(meta_to_carbon)'
	substructure = '[#6:1]~[a:2]1[a;!$(a(a)(a)[#6,#7,#8]):3][cH1:4][a;!$(a(a)(a)[#6,#7,#8]):5][a:6][a;!$(a(a)(a)[#6,#7,#8]):7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_6(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(ortho_to_nitrogen)'
	substructure = '[#7:1]~[a:2]1[cH1:3][a;!$(a(a)(a)[#6,#7,#8]):4][a:5][a;!$(a(a)(a)[#6,#7,#8]):6][a:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_7(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(ortho_to_oxygen)'
	substructure = '[#8:1]~[a:2]1[cH1:3][a;!$(a(a)(a)[#6,#7,#8]):4][a:5][a;!$(a(a)(a)[#6,#7,#8]):6][a:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_8(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(ortho_to_2_substituents)'
	substructure = '[#6,#7,#8:1]~[a:2]1[cH1:3][a;$(a(a)(a)[#6,#7,#8]):4][a:5][a;!$(a(a)(a)[#6,#7,#8]):6][a:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_9(molecule: Chem.rdchem.Mol):
	# 'aromatic_hydroxylation_(sulfur_containing_5ring)'
	substructure = '[cH1;$(c1saaa1):2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_10(molecule: Chem.rdchem.Mol):
	# 'aromatic_oxidation_(nitrogen_containing_5ring)'
	substructure = '[nH0:1][cH1;$(c1naan1):2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_ aromatic hydroxylation_11(molecule: Chem.rdchem.Mol):
	# 'aromatic_dehydroxylation'
	substructure = '[c;$(cc[OH1]):1][OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_carboxylation_1(molecule: Chem.rdchem.Mol):
	# 'carboxylation_(primary_carbon_next_to_quart_carbon)'
	substructure = '[C;X4;H0;$(C[!C]):1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_carboxylation_2(molecule: Chem.rdchem.Mol):
	# 'carboxylation_(primary_carbon_next_to_tert_carbon)'
	substructure = '[CH1;$(C(-[#6])(-[#6])-[CH3]):1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_carboxylation_3(molecule: Chem.rdchem.Mol):
	# 'carboxylation_(primary_carbon_next_to_sec_carbon)'
	substructure = '[#6:1][CH2:2][CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_carboxylation_4(molecule: Chem.rdchem.Mol):
	# 'carboxylation_(primary_carbon_next_to_SP2)'
	substructure = '[C;$(C=*),$(C#*):1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_carboxylation_5(molecule: Chem.rdchem.Mol):
	# 'carboxylation_(benzylic_CH3)'
	substructure = '[c:1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_1(molecule: Chem.rdchem.Mol):
	# 'general_all_aliph_hydr'
	substructure = '[C;X4;!H0;!$(Cc):1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_2(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(primary_carbon_next_to_quart_carbon)'
	substructure = '[C;X4;H0;$(C[!C]):1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_3(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(primary_carbon_next_to_tert_carbon)'
	substructure = '[CH1;$(C(-[#6])(-[#6])-[CH3]):1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_4(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(primary_carbon_next_to_sec_carbon)'
	substructure = '[#6:1][CH2:2][CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_5(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(primary_carbon_next_to_SP2_or_SP1)'
	substructure = '[C;$(C=*),$(C#*):1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_6(molecule: Chem.rdchem.Mol):
	# 'unspecific_secondary_aliphatic_carbon_hydroxylation'
	substructure = '[CX4:1][CH2:2][CX4:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_7(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(sec_carbon,next_to_CH3)'
	substructure = '[CX4:1][CH2:2][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_8(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(sec_carbon_in_a_ringA)'
	substructure = '[CX4;H2:1][CH2;R:2][CX4;H2:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_9(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(sec_carbon_in_a_ringB)'
	substructure = '[CX4;H2:1][CH2;R:2][CX4;!H2:3][*;$([CH3]),!#6:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_10(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(sec_carbon_next_to_SP2,not_in_a_ring)'
	substructure = '[CX4:1][CH2;!R:2][*;!c;$(*=*):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_11(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(sec_carbon_next_to_SP2,in_a_ring)'
	substructure = '[CX4:1][CH2;R:2][*;!c;$(*=*),$([#7]):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_12(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(sec_carbon_both_sides_next_to_SP2,in_a_ring)'
	substructure = '[*;!c;$(*=*):1][CH2;R:2][*;!c;$(*=*):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_13(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(tert_carbon_next_to_SP2)'
	substructure = '[C:1][CH1;X4:2]([C;!$([CH3]):3])[N,C&$([C]=*):4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aliphatic_hydroxylation_14(molecule: Chem.rdchem.Mol):
	# 'aliphatic_hydroxylation_(tert_carbon_linked_to_two_CH3_groups)'
	substructure = '[CH3][CH1;X4;!$(Cc):1][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_benzylic_hydroxylation_1(molecule: Chem.rdchem.Mol):
	# 'benzylic_hydroxylation_(c-CH3)'
	substructure = '[c:1][CH3:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_benzylic_hydroxylation_2(molecule: Chem.rdchem.Mol):
	# 'benzylic_hydroxylation_(c-CH2-CH3)'
	substructure = '[c:1][CH2:2][CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_benzylic_hydroxylation_3(molecule: Chem.rdchem.Mol):
	# 'benzylic_hydroxylation_(c-CH2-CR)'
	substructure = '[c:1][CH2:2][#6;!$([CH3]):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_benzylic_hydroxylation_4(molecule: Chem.rdchem.Mol):
	# 'benzylic_hydroxylation_(c-CH2-N)'
	substructure = '[c:1][CH2:2][NH0:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_benzylic_hydroxylation_5(molecule: Chem.rdchem.Mol):
	# 'benzylic_hydroxylation_(c-CH1-CH3)'
	substructure = '[c:1][CH1;X4;!$(C[O,N]):2][CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_benzylic_hydroxylation_6(molecule: Chem.rdchem.Mol):
	# 'benzylic_hydroxylation_(c-CH1-CR)'
	substructure = '[c:1][CH1;X4;!$(C[O,N]):2][#6;c,$(C=*):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_1(molecule: Chem.rdchem.Mol):
	# 'carbonyl_reduction_(aliphatic)'
	substructure = '[C;X4:1][C:2](=[O:3])[C;X4:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_2(molecule: Chem.rdchem.Mol):
	# 'carbonyl_reduction_(next_to_SP2_carbon)'
	substructure = '[C;X3:1][C:2](=[O:3])[C;X4:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_3(molecule: Chem.rdchem.Mol):
	# 'carbonyl_reduction_(next_to_aromatic_carbon)'
	substructure = '[c:1][C:2](=[O:3])[C;X4:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_4(molecule: Chem.rdchem.Mol):
	# 'carbonyl_reduction_(both_sides_next_to_aromatic_carbon)'
	substructure = '[c:1][C:2](=[O:3])[c:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_5(molecule: Chem.rdchem.Mol):
	# 'aldehyde_reduction_(aliphatic)'
	substructure = '[C:1][CH1:2]=[O:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_6(molecule: Chem.rdchem.Mol):
	# 'aldehyde_reduction_(aromatic)'
	substructure = '[c:1][CH1:2]=[O:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_7(molecule: Chem.rdchem.Mol):
	# 'double_bond_reduction'
	substructure = '[C;$(C[OH1]),$(C=O):1][C:2]=[C;!$(Cc):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_8(molecule: Chem.rdchem.Mol):
	# 'double_bond_reduction_(aromatic)'
	substructure = '[c;$(c=O):1][c:2][cH1;$(co),$(cn):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_reduction_9(molecule: Chem.rdchem.Mol):
	# 'double_bond_reduction_(benzylic)'
	substructure = '[C;$(C[OH1]),$(C=O):1][C:2]=[C;$(Cc):3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aldehyde_oxidation_1(molecule: Chem.rdchem.Mol):
	# 'aldehyde_oxidation_(aliphatic)'
	substructure = '[C:1][CH1:2]=[O:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_aldehyde_oxidation_2(molecule: Chem.rdchem.Mol):
	# 'aldehyde_oxidation_(aromatic)'
	substructure = '[c:1][CH1:2]=[O:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_O_deacetylation_1(molecule: Chem.rdchem.Mol):
	# 'O-deacetylation'
	substructure = '[#6:1][O:2]C(=O)[CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_deacetylation_1(molecule: Chem.rdchem.Mol):
	# 'N-deacetylation'
	substructure = '[N:2]C(=O)[CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_decarboxylation_1(molecule: Chem.rdchem.Mol):
	# 'decarboxylation'
	substructure = '[*;!C:1]~[#6:2]C(=O)[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_decarboxylation_2(molecule: Chem.rdchem.Mol):
	# 'oxidative_decarboxylation'
	substructure = '[O:1]=[C:2][C:3](=O)[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_decarboxylation_3(molecule: Chem.rdchem.Mol):
	# 'beta-oxidation'
	substructure = '[CH2:1][CH2]C(=O)[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_1(molecule: Chem.rdchem.Mol):
	# 'general'
	substructure = '[C:1][C:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_2(molecule: Chem.rdchem.Mol):
	# 'dehydrogenation_(alpha,beta_to_SP2_both_sides)'
	substructure = '[*;$([#6&X3]),$([#7]~[#6X3]):1][CX4;H1&!$(C-[!#6]),H2:2][CX4;H2:3][*;$([#6&X3]),$([#7]~[#6X3]):4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_3(molecule: Chem.rdchem.Mol):
	# 'dehydrogenation_(alpha,beta_to_SP2)'
	substructure = '[*;$([#6&X3]),$([#7]~[#6X3]):1][CX4;H1&!$(C-[!#6]),H2:2][CX4;H2:3][C;H2,H3:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_4(molecule: Chem.rdchem.Mol):
	# 'dehydrogenation_(CH1-CH3->C=CH2)'
	substructure = '[#6X3:1][CH1&!$(C-[!#6]):2][CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_5(molecule: Chem.rdchem.Mol):
	# 'dehydrogenation_(CH2-CH3->C=CH2)'
	substructure = '[#6X3:1][CH2:2][CH3:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_6(molecule: Chem.rdchem.Mol):
	# 'dehydrogenation_(amine)'
	substructure = '[N,c:1][C;X4;H1:2]-[N;X3;H1:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydrogenation_7(molecule: Chem.rdchem.Mol):
	# 'dehydrogenation_(aromatization_of_1,4-dihydropyridine)'
	substructure = '[c:1][#6:2]1[#6:3]=[#6:4][NH1:5][#6:6]=[#6:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydration_1(molecule: Chem.rdchem.Mol):
	# 'dehydration_next_to_SP2_both_sides'
	substructure = '[CX4@!H0;$(C[*;#6&X3,$([#7]~[#6X3])]):1]-[CX4@;$(C[*;#6&X3,$([#7]~[#6X3])]):2]([OH1])'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydration_2(molecule: Chem.rdchem.Mol):
	# 'dehydration_next_to_SP2_a'
	substructure = '[CX4@!H0;!$(C[*;#6&X3,$([#7]~[#6X3])]):1]-[CX4@;$(C[*;#6&X3,$([#7]~[#6X3])]):2]([OH1])'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehydration_3(molecule: Chem.rdchem.Mol):
	# 'dehydration_next_to_SP2_b'
	substructure = '[CX4@!H0;$(C[*;#6&X3,$([#7]~[#6X3])]):1]-[CX4@;!$(C[*;#6&X3,$([#7]~[#6X3])]):2]([OH1])'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_primary_alcohol_oxidatyion_to_carboxyl_1(molecule: Chem.rdchem.Mol):
	# 'primary_alcohol_oxidation_(benzylic)'
	substructure = '[c:1][CH2:2][OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_primary_alcohol_oxidatyion_to_carboxyl_2(molecule: Chem.rdchem.Mol):
	# 'primary_alcohol_oxidation_(aliphatic)'
	substructure = '[C:1][CH2:2][OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_secondary_alcohol_oxidation_to_carbonyl_1(molecule: Chem.rdchem.Mol):
	# 'secondary_alcohol_oxidation_(aliphatic)'
	substructure = '[C;!$(C[OH1]):1][CH1:2]([C;!$(C[OH1]):3])-[OH1:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_secondary_alcohol_oxidation_to_carbonyl_2(molecule: Chem.rdchem.Mol):
	# 'secondary_alcohol_oxidation_(benzylic)'
	substructure = '[c:1][CH1:2]([C:3])-[OH1:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_1(molecule: Chem.rdchem.Mol):
	# 'sulfoxide_oxidation_(c-S-C)'
	substructure = '[c:1][S;X3:2](=[O:3])[C:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_2(molecule: Chem.rdchem.Mol):
	# 'sulfoxide_oxidation_(C-S-C)'
	substructure = '[C:1][S;X3:2](=[O:3])[C:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_3(molecule: Chem.rdchem.Mol):
	# 'sulfoxide_oxidation_(c-S-c)'
	substructure = '[c:1][S;X3:2](=[O:3])[c:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_4(molecule: Chem.rdchem.Mol):
	# 'sulfide_oxidation_(c-S-C)'
	substructure = '[c:1][S;X2:2][C:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_5(molecule: Chem.rdchem.Mol):
	# 'sulfide_oxidation_(C-S-C)'
	substructure = '[C:1][S;X2:2][C:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_6(molecule: Chem.rdchem.Mol):
	# 'sulfide_oxidation_(c-S-c)'
	substructure = '#[c:1][#16;X2:2][c:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_7(molecule: Chem.rdchem.Mol):
	# 'sulfide_oxidation_(c-S-c)'
	substructure = '[c:1][S;X2:2][c:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_8(molecule: Chem.rdchem.Mol):
	# 'thiophene_oxidation'
	substructure = '[sr5:1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_S_oxidation_9(molecule: Chem.rdchem.Mol):
	# 'sulfoxide_reduction'
	substructure = '[S;X3;$(S([#6])[#6]):1]=O'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_epoxide_hydrolysis_1(molecule: Chem.rdchem.Mol):
	# 'epoxide_hydrolysis'
	substructure = '[C:1]1O[C:2]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_epoxide_hydrolysis_2(molecule: Chem.rdchem.Mol):
	# 'epoxidation'
	substructure = '[C:1]=[C:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_oxidative_deamination_1(molecule: Chem.rdchem.Mol):
	# 'oxidative_deamination_(amidine)'
	substructure = '[#6:1][N:2]=;@[C:3]([#6:4])[N:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_oxidative_deamination_2(molecule: Chem.rdchem.Mol):
	# 'oxidative_deamination_(aromatic)'
	substructure = '[nX2:1][c:2][N:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_oxidative_deamination_3(molecule: Chem.rdchem.Mol):
	# 'oxidative_deamination_(on_primary_carbon)'
	substructure = '[C:1][CH2:2][NH2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_oxidative_deamination_4(molecule: Chem.rdchem.Mol):
	# 'oxidative_deamination_(on_secondary_carbon)'
	substructure = '[C:1][CH1:2]([C:3])[NH2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_nitro_1(molecule: Chem.rdchem.Mol):
	# 'nitro_to_aniline'
	substructure = '[c:1][N+](=O)[O-]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_nitro_2(molecule: Chem.rdchem.Mol):
	# 'aniline_to_nitro'
	substructure = '[c;$(c1[cH1][cH1][c]([*;!#1])[cH1][cH1]1):1][NH2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_nitro_3(molecule: Chem.rdchem.Mol):
	# 'nitro_to_nitroso'
	substructure = '[c:1][N+:2](=[O])[O-]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehalogenation_1(molecule: Chem.rdchem.Mol):
	# 'haloacid_hydrolysis'
	substructure = '[#6:1][C:2](=[O:3])[*;F,Cl,Br,I]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehalogenation_2(molecule: Chem.rdchem.Mol):
	# 'oxidative_dehalogenation'
	substructure = '[C:1]([OH1:2])[*;Cl,Br,I]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehalogenation_3(molecule: Chem.rdchem.Mol):
	# 'aliphatic_dehalogenation'
	substructure = '[CX4;H1,H2:1][Cl,Br,I]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_dehalogenation_4(molecule: Chem.rdchem.Mol):
	# 'aromatic_dechlorination'
	substructure = '[c;$(c1ccc([#7])cc1):1][Cl]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_1a(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-5bonds-carboxyl)'
	substructure = '[OH1][C:2]!@[*:3]~!@[*:4][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_1b(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-5bonds-carboxyl)'
	substructure = '[OH1][C:2]@[*:3]~!@[*:4][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_1c(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-5bonds-carboxyl)'
	substructure = '[OH1][C:2]!@[*:3]~@[*:4][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_2a(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-5bonds-carboxyl)2'
	substructure = '[NH1;!$(NC=O):1][#6:2]~!@[*:3]~!@[*:4]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_2b(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-5bonds-carboxyl)2'
	substructure = '[NH1;!$(NC=O):1][#6:2]~[*:3]~!@[*:4]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_2c(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-5bonds-carboxyl)2'
	substructure = '[NH1;!$(NC=O):1][#6:2]~!@[*:3]~[*:4]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_3a(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-6bonds-carboxyl)'
	substructure = '[OH1][C:2]!@[*:3]~!@[*:4]~!@[*:5][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_3b(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-6bonds-carboxyl)'
	substructure = '[OH1][C:2]@[*:3]~!@[*:4]~!@[*:5][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_3c(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-6bonds-carboxyl)'
	substructure = '[OH1][C:2]!@[*:3]~@[*:4]~!@[*:5][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_3d(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(hydroxyl-6bonds-carboxyl)'
	substructure = '[OH1][C:2]!@[*:3]~!@[*:4]~@[*:5][C;!$(CC1OCC(O)C(O)C1O)](=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_5a(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-6bonds-carboxyl)'
	substructure = '[NH1;!$(NC=O):1][#6:2]~!@[*:3]~!@[*:4]~!@[*:5]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_5b(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-6bonds-carboxyl)'
	substructure = '[NH1;!$(NC=O):1][#6:2]~@[*:3]~!@[*:4]~!@[*:5]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_5c(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-6bonds-carboxyl)'
	substructure = '[NH1;!$(NC=O):1][#6:2]~!@[*:3]~@[*:4]~!@[*:5]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_5d(molecule: Chem.rdchem.Mol):
	# 'ring_closure_(NH1-6bonds-carboxyl)'
	substructure = '[NH1;!$(NC=O):1][#6:2]~!@[*:3]~!@[*:4]~@[*:5]C(=O)-[OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_6a(molecule: Chem.rdchem.Mol):
	# 'hydroxyl-amide_5ring_closure'
	substructure = '[OH1:1][C:2][A:3][A:4][C:5](=[O:6])[N:7]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_6b(molecule: Chem.rdchem.Mol):
	# 'hydroxyl-amide_6ring_closure'
	substructure = '[OH1:1][C:2][A:3][A:4][A:5][C:6](=[O:7])-[N:8]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_6c(molecule: Chem.rdchem.Mol):
	# 'hydroxyl-amide_5ring_rearr'
	substructure = '[OH1:1][C:2][A:3][N:4][C:5](=[O:6])'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_condensation_6d(molecule: Chem.rdchem.Mol):
	# 'hydroxyl-amide_6ring_rearr'
	substructure = '[OH1:1][C:2][A:3][A:4][N:5][C:6](=[O:7])'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_1(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(methoxyester)'
	substructure = '[C;$(C=O):1][O:2][CH3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_2(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(ester)'
	substructure = '[C$(C[#6!H3]):2](=[O:3])O[#6!H3:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_3(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(primary_amide)'
	substructure = '[C$(C[#6!H3]):2](=[O:3])[NH2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_4(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(secondary_amide)'
	substructure = '[C$(C[#6!H3]):2](=[O:3])[NH1:4][#6:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_5(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(tertiary_amide)'
	substructure = '[C$(C[#6!H3]):2](=[O:3])[#7:4]([#6:5])[#6:6]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_6(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(heteroatom_bonded_amide)'
	substructure = '[C$(C[#6!H3]):2](=[O:3])[N:4][*;!#6:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_7(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(urea_or_carbonate)'
	substructure = '[#7,#8:1][C:2](=[O:3])[#7,#8:4][*:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_8(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(X=X-X_exclude_phosphate)'
	substructure = '[*:5][*;!#6;!$(S(=O)(=O)N);!$(P(O)(O)(O)=O):1](=[*;!#6:2])[N,O:3][*:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_9(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(CNC(OH)R)'
	substructure = '[#6:1][N:2][CH1]([OH1])[*:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_hydrolysis_10(molecule: Chem.rdchem.Mol):
	# 'hydrolysis_(N-substituted-pyridine)'
	substructure = '[n:2][c:3]!@[N;$(N(C)(C)c),$(NS(=O)=O):5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_oxidation_1(molecule: Chem.rdchem.Mol):
	# 'N-oxidation_(tertiary_N)'
	substructure = '[C;X4;!H3;!$(C(N)[!#6;!#1]):1][N;X3:2]([C;X4;!H3;!$(C(N)[!#6;!#1]):3])[C;X4;!H3;!$(C(N)[!#6;!#1]):4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_oxidation_2(molecule: Chem.rdchem.Mol):
	# 'N-oxidation_(tertiary_NCH3)'
	substructure = '[C;X4;!H3;!$(C(N)[!#6;!#1]):1][N;X3:2]([CH3:3])[C;X4;!H3;!$(C(N)[!#6;!#1]):4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_oxidation_3(molecule: Chem.rdchem.Mol):
	# 'N-oxidation_(RN(CH3)2)'
	substructure = '[C;X4;!$(C(N)[!#6;!#1]):1][N;X3:2]([CH3:3])[CH3:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_oxidation_4(molecule: Chem.rdchem.Mol):
	# 'N-oxidation_(-N=)'
	substructure = '[#6:1]~[#7;X2;R:2]~[#6:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_N_oxidation_5(molecule: Chem.rdchem.Mol):
	# 'N-oxidation_(aniline)'
	substructure = '[c:1][NH2:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_acetyl_shift_1(molecule: Chem.rdchem.Mol):
	# 'acetyl_shift'
	substructure = '[#6:1][C:2](=O)O[C:5][C:6][OH1]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_tautomerisation_1(molecule: Chem.rdchem.Mol):
	# 'tautomerisation_(keto->enol)'
	substructure = '[c:1][C:2](=[O:3])[CH2:4][#6:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_1(molecule: Chem.rdchem.Mol):
	# 'vinyl_oxidation'
	substructure = '[#6:3][CH1:1]=[CH2:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_2(molecule: Chem.rdchem.Mol):
	# 'isopropenyl_oxidation'
	substructure = '[#6:3][C:1]([CH3:4])=[CH2:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_3(molecule: Chem.rdchem.Mol):
	# 'oxidation_(amine_in_a_ring)'
	substructure = '[CH2:1][CH2;R:2][N:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_4(molecule: Chem.rdchem.Mol):
	# 'imine_hydrolysis'
	substructure = '[#6:1][C:2]([#6:3])=[N;!$(N-N):4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_5(molecule: Chem.rdchem.Mol):
	# 'hydrazone_hydrolysis'
	substructure = '[#6:2]=[N:4]-[N:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_6(molecule: Chem.rdchem.Mol):
	# 'diazene_cleavage'
	substructure = '[c:1][N:2]=[N:3][c:4]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_7(molecule: Chem.rdchem.Mol):
	# 'azide_cleavage'
	substructure = '[*:1][N:2]=[N+]=[N-]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_8(molecule: Chem.rdchem.Mol):
	# 'aromatic_oxidation'
	substructure = '[#6:1][c:2]1[cH1:3][cH1:4][cH1:5][cH1:6][cH1:7]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_9(molecule: Chem.rdchem.Mol):
	# 'phosphine_sulphide_hydrolysis'
	substructure = '[P:1]=[S]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_10(molecule: Chem.rdchem.Mol):
	# 'xanthine_oxidation'
	substructure = '[N:1][CH1:2]=[N:3]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_11(molecule: Chem.rdchem.Mol):
	# 'oxidation_to_quinone'
	substructure = '[#7,O;H1:1][#6:2]:1:[#6:3]:[#6:4]:[#6:5](:[#6:6]:[#6:7]:1)[#7,O;H1:8]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_12(molecule: Chem.rdchem.Mol):
	# 'cyclic_hemiacetal_ring_opening'
	substructure = '[#6:1][O:2]@[CH1:3]([OH1:4])[*:5]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_13(molecule: Chem.rdchem.Mol):
	# 'oxidation_(C=N)'
	substructure = '[NX2:1]=[CH1:2]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_14(molecule: Chem.rdchem.Mol):
	# 'deiodonidation'
	substructure = '[#6X3:1][I]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_special_rules_15(molecule: Chem.rdchem.Mol):
	# 'nitrile_to_amide'
	substructure = '[C:1]#N'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_steroids_1(molecule: Chem.rdchem.Mol):
	# 'steroid_d5d4'
	substructure = '[C;$(C~1~C~C~C~C~2~C~C~C~3~C~4~C~C~C~C~4~C~C~C~3~C~2~1):1]1[C:2][C:3](=[O:30])[C:4][C:5]=[C:6]1'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)

def REACTION_steroids_2(molecule: Chem.rdchem.Mol):
	# 'steroid_17hydroxy_to_keto'
	substructure = '[C;$(C~1~C~2~C~C~C~3~C~4~C~C~C~C~C~4~C~C~C~3~C~2~C~C~1):17]([OH1:30])!@[C:31]'
	substructure = Chem.MolFromSmarts(substructure)
	indices = molecule.GetSubstructMatches(substructure)
	return len(indices)
