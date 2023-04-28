import csv
import numpy as np
import matplotlib.pyplot as plt
from urllib.request import urlopen
import pandas as pd
import openpyxl
from sbtab import SBtab

import equilibrator_api
import math


from equilibrator_api import ComponentContribution, Q_, R, default_T
import equilibrator_pathway
from equilibrator_pathway import Pathway
from equilibrator_pathway.thermo_models import ThermodynamicModel
from equilibrator_pathway.mdf_solution import PathwayMdfSolution
from equilibrator_pathway.cost_function import EnzymeCostFunction
from equilibrator_pathway.ecm_model import EnzymeCostModel
from equilibrator_api.model.bounds import Bounds
from equilibrator_api.phased_reaction import PhasedReaction


from equilibrator_pathway.cost_function import EnzymeCostFunction
from optlang.glpk_interface import Constraint, Model, Objective, Variable
from scipy.optimize import minimize
import logging
from IPython.display import display, HTML

"""
FUNCTIONS FOR BASIC READING OF CSV FILES
"""
cc = ComponentContribution()
RT = R * default_T

#Class which packages all objects obtained from csv metabolite file (with obtain_compound_settings)
class CompoundSettings(object):
	def __init__(self, sd, nd, ad, bs):#,conserved_moieties = []):
		self.acr2kegg_dict = sd #Acronym to keggID dictionary
		self.cmpstrs2acr_dict = nd #Compound in string form to acronym dictionary
		self.cmpd2acr_dict = ad #Compound to acronym dictionary
		if bs is None:
			bs = Bounds.get_default_bounds(cc).copy()       # old - GetDefaultBounds().Copy()
		self.bounds = bs
		self.acr2cmpd_dict = {value : key for (key, value) in ad.items()} #Acronym to compound dictionary
		#self.conserved_moieties = [ConservedMoiety(acr_stoi,total_conc,self) for (acr_stoi,total_conc) in conserved_moieties]

	#Copy CompoundSettings object
	def Copy(self):
		return CompoundSettings(self.acr2kegg_dict.copy(),self.cmpstrs2acr_dict.copy(),self.cmpd2acr_dict.copy(),self.bounds.copy())#,conserved_moieties=self.conserved_moieties) # old-bounds.Copy

	#Change bounds of CompoundSettings in-place
	def change_bounds(self,bound_sets):
		if type(bound_sets) is not list: bound_sets = [ bound_sets ]
		for (ca, lb, ub) in bound_sets:
			compound = self.acr2cmpd_dict.get(ca)
			if compound is None:
				print(arg + " not a valid acronym!")
			else:
				self.bounds.set_bounds(compound, lb, ub) #old-SetBounds
	
	#Obtain lower bound, upper bound or bound tuple from CompoundSettings object
	def get_bounds(self,cmpd_acronym, bound_type = 'both'):
		if bound_type == 'both': # old-bound_type is 'both':
			return self.bounds.get_bound_tuple(self.acr2cmpd_dict[cmpd_acronym])
		if bound_type == 'lower':
			return self.bounds.get_lower_bounds(self.acr2cmpd_dict[cmpd_acronym])
		if bound_type == 'upper':
			return self.bounds.get_upper_bound(self.acr2cmpd_dict[cmpd_acronym])


#Takes metabolite reference file to produce dictionaries for converting between compound IDs,
#as well as setting up boundaries for each 

def obtain_compound_settings(filename,custom_bounds = True):
	with open(filename + '.csv','r') as csv_file:
		df = pd.read_csv(csv_file)
		csv_file.close()
	acronyms = df['Metabolite Acronym']
	#print(acronyms)
	keggIDs = [k.strip() for k in df['keggID']]
	compounds = list(map(cc.get_compound,keggIDs))  #equilibrator_api.ccache.get_compound
	
	if custom_bounds:
		l_bds = [Q_(x, 'M') for x in map(float,df['Lower Bound (M)'])]
		u_bds = [Q_(x, 'M') for x in map(float,df['Upper Bound (M)'])]
		bds = equilibrator_pathway.Bounds()
		for i in range(len(compounds)):
			bds.set_bounds(compounds[i],l_bds[i],u_bds[i])    #old-SetBounds
	else:
		bds = None
	#print(bds.GetBounds(compounds))
	#dict(zip(keys, values))         
	acr2kegg_dict = dict(zip(acronyms,keggIDs))
	cmp_strs = [c.__repr__() for c in compounds]
	cmpstrs2acr_dict = dict(zip([cs.replace(', ',',',1).replace('\n','',1) for cs in cmp_strs],acronyms))
	cmpd2acr_dict = dict(zip(compounds,acronyms))
	
	return CompoundSettings(acr2kegg_dict,cmpstrs2acr_dict,cmpd2acr_dict,bds)


# my-Ergänzung partly*
#Takes pathway csv file, compound settings and physiological context to translate reactions and fluxes for the api,
#also naming reactions with the given names, returning the pathway object. Also, if standard dG is provided, it is implemented in the pathway, unless it is stated otherwise
def create_pathway(filename, sheetname,cs,pc,absolute_flux = Q_(1.0, 'M/s'),custom_dGs = True,x = 0):
	(reactions, fluxes,std_dGs) = custom_parse_reactions(filename,sheetname,cs,absolute_flux = absolute_flux,custom_dGs = custom_dGs)
	pathway = equilibrator_pathway.ThermodynamicModel(reactions, fluxes,standard_dg_primes=std_dGs,bounds=cs.bounds,config_dict = pc)   # equilibrator_pathway.Pathway
		
	pathway.dg_confidence = x

	pathway.set_compound_ids(CallAcronyms(cs.cmpd2acr_dict))                                                                                
	return pathway

#Converts data from pathway csv file into a list of Reaction objects, an array of relative fluxes and custom dG values, if defined so
def custom_parse_reactions(filename,sheetname,cs,absolute_flux = Q_(1.0, 'M/s'),custom_dGs = True):
	with open(filename + '.xlsx','r', encoding='UTF-8') as excel_file:
		df = pd.read_excel(excel_file.name,sheet_name=sheetname)
		#print(df)
		excel_file.close()
	rxns = df.get('Reaction Formula')
	if rxns is None:
		print('Reaction Formula column not found!')
	#fluxes = np.array(df['RelativeFlux']) * Q_("dimensionless")
	#flux_strings = df.get('AbsoluteFlux(M/s)')
	#if flux_strings is None:
	flux_strings = df.get('Relative Flux')
	if flux_strings is None:
		print('Relative Flux column not found!')
	fluxes = np.array(flux_strings, ndmin=2, dtype=float).T
	fluxes *= absolute_flux
	rids = df.get('Reaction Name')
	if rids is None:
		print('Reaction Name column not found!')
	new_rxns = []
	for i in range(len(rxns)):
		new_rxn = ""
		#print(rxns[i])
		args = rxns[i].split(' ')
		for arg in args:
			#print( arg )
			ANumber = all(map(str.isdigit, arg.replace(".", "", 1) ))
			if(arg != '+' and arg != '<=>' and arg != '=' and not ANumber):
				keggID = cs.acr2kegg_dict.get(arg)
				if keggID == None:
					print(arg + " not a valid acronym!")
				else:
					arg = keggID
			new_rxn += arg + ' '
		new_rxn = PhasedReaction.parse_formula(cc.get_compound, new_rxn[:-1],rids[i])
		if not new_rxn.is_balanced():
			print('%s is not balanced' % new_rxn.rid)
			print( new_rxn._get_reaction_atom_bag() )
		new_rxns.append(new_rxn)
	reactions = new_rxns
	std_dGs = None
	if custom_dGs:
		std_dGs = df.get('Standard dG (kJ/mol)')
		if (std_dGs is not None):
			std_dGs = np.array(std_dGs, ndmin=2, dtype=float).T.flatten()
			std_dGs *= Q_(1.0,'kJ/mol')

	return (reactions,fluxes,std_dGs)

#Callable class of metabolite acronyms, specific for use as argument in pathway.set_compound_ids function
class CallAcronyms:
	def __init__(self, cmpd2acr_dict):
		self.cmpd2acr_dict = cmpd2acr_dict
	def __call__(self,compound):
		return self.cmpd2acr_dict[compound]

#Makes reaction string with compound object ids, into reaction string with specified acronyms
def translate_reaction_str(reaction,cmpstrs2acr_dict):
	reaction = reaction.replace(', ',',').replace('\n','')
	reaction = reaction.split()
	transl_rxn = ""
	for elem in reaction:
		if elem.startswith('Compound') and elem in cmpstrs2acr_dict:
			transl_rxn += cmpstrs2acr_dict[elem] + " "
		else:
			transl_rxn += elem + " "
	return transl_rxn[:-1]

#Makes reaction object into reaction string with specified acronyms
def translate_reaction(reaction,cmpstrs2acr_dict):
	reaction = reaction.__str__()
	return translate_reaction_str(reaction,cmpstrs2acr_dict)


"""
FUNCTIONS AND STRUCTURES TO EDIT ANALYSIS
"""

#Returns a CompoundSettings object with updated bounds.
#The bound sets should be a list of tuples or a single tuple like: (str metaboliteAcronym, Quantity lowerBound, Quantity higherBound)
def change_bounds(cs : CompoundSettings,bound_sets):
	new_cs = cs.Copy()
	new_cs.change_bounds(bound_sets)
	return new_cs

#Return a CompoundSettings object with fixed concentrations
#Same as change_bounds, but input need only to be a duple like (str metaboliteAcronym, Quantity both_bounds)
def fix_concentrations(cs, conc_sets):
	if type(conc_sets) is not list: conc_sets = [ conc_sets ]
	return change_bounds(cs, [(acr, conc, conc) for (acr,conc) in conc_sets])

#Changes bounds of defined compounds to default eQuilibrator metabolite bounds, 1 uM <= C <= 10 mM
def free_concentrations(cs, acronyms):
	if type(acronyms) is not list: acronyms = [ acronyms ]
	return change_bounds(cs,[(acr,Q_(1e-6,'M'),Q_(1e-2,'M')) for acr in acronyms])


"""
FUNCTIONS FOR MDF ANALYSIS
"""

#Full pre-processing, configuring, calculation and plotting of MDF
def MDF_analysis(pathway_file_name, pathway_sheet_name, cmpd_settings, physiological_context, custom_dGs = True, print_results = True, y = 0):
	#compound = cmpd_settings.acr2cmpd_dict.get("AcCoA")
	
	pathway = create_pathway(pathway_file_name, pathway_sheet_name,cmpd_settings,physiological_context,custom_dGs=custom_dGs, x = y)
	
	return MDF_from_pathway(pathway, cmpd_settings, print_results = print_results)
	#print(list(cmpstrs2acr_dict.keys())[list(cmpstrs2acr_dict.values()).index('F6P')])

#Performs MDF analysis from pathway object, can be used with create pathway to complete
def MDF_from_pathway(pathway, cmpd_settings, print_results = True):
	mdf_result = pathway.mdf_analysis()         # pathway.calc_mdf()
	mdf_result.reaction_df["reaction_formula"] = [translate_reaction(rxn,cmpd_settings.cmpstrs2acr_dict) for rxn in mdf_result.reaction_df["reaction_formula"]]

	if not pathway.net_reaction.is_balanced():
		print("Net reaction is not balanced!")
	if print_results:
		print(translate_reaction(pathway.net_reaction,cmpd_settings.cmpstrs2acr_dict))
		print_MDF_results(mdf_result)
	return mdf_result


# Displays a full report of MDF results with compound concentration plot and and cumulative dG plot of reactions.
#Cmpd_settings argument should be used only in the case of the result being printed originates from multi-result functions, such as ratio_sweep,
#as it is used to translate reaction formulae with the user defined compound names.
def print_MDF_results(mdf_result, cmpd_settings = None):
	if cmpd_settings is not None:
		mdf_result.reaction_df["reaction_formula"] = [translate_reaction(rxn,cmpd_settings.cmpstrs2acr_dict) for rxn in mdf_result.reaction_df["reaction_formula"]]
		#mdf_result.compound_df["compound"] = [cmpd_settings.cmpstrs2acr_dict.get(cmpd) for cmpd in mdf_result.compound_df]
	display(f"MDF = {mdf_result.score:.2f} kJ/mol")  # old- f"MDF = {mdf_result.mdf:.3f}"
	#display(mdf_result.reaction_df)
	#display(mdf_result.compound_df)
	#fig1 = mdf_result.plot_concentrations();      # old-mdf_result.compound_plot
	#fig2 = mdf_result.plot_driving_forces();


# my_Ergänzung
def print_MDF_results_individual(mdf_result):
	
	MDF  = mdf_result.score
  	
	res1 = mdf_result.reaction_df
	g1   = mdf_result.reaction_df["original_standard_dg_prime"]
	g2   = mdf_result.reaction_df["standard_dg_prime"]
	g3   = mdf_result.reaction_df["physiological_dg_prime"]
	g4   = mdf_result.reaction_df["optimized_dg_prime"]
	flx  = mdf_result.reaction_df["flux"]
	shp_r= mdf_result.reaction_df["shadow_price"]
	r_id = mdf_result.reaction_df["reaction_id"]

	res2 = mdf_result.compound_df
	met  = mdf_result.compound_df["compound_id"]
	bl   = mdf_result.compound_df["lower_bound"]
	bu   = mdf_result.compound_df["upper_bound"]
	c    = mdf_result.compound_df["concentration"]
	shp_c= mdf_result.compound_df["shadow_price"]
	
	return MDF,res1,g1,g2,g3,g4,flx,shp_r,r_id,res2,met,bl,bu,c,shp_c
	

#Obtains the shadow price of a single compound
def cmpd_shadow_price(mdf_result,cmpd_acronym):
	shadow_price = mdf_result.compound_df.loc[mdf_result.compound_df["compound"] == cmpd_acronym]["shadow_price"]
	shadow_price = shadow_price.values[0] #Get the value from the returned series
	return shadow_price

#Calculates the MDF across a range of concentrations (in Molar) of a specific metabolite in the pathway
def MDF_conc_sweep(compound_acronym, conc_range,pathway_file_name, compound_settings, physiological_context, custom_dGs = True, y = True):
	cs = compound_settings
	
	if y == True:
		pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x = True)
	else:
		pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x = y)
	
	mdfs = []
	compound = cs.acr2cmpd_dict.get(compound_acronym)
	if compound is None:
		print(compound_acronym + " not a valid acronym!")
	for i in range(conc_range.size):
		pathway._bounds.set_bounds(compound, Q_(conc_range[i], "M"), Q_(conc_range[i], "M"))	
		mdf_result = pathway.mdf_analysis()
		mdf_Q = mdf_result.score
		mdfs.append(mdf_Q)
	#mdfs = np.array([mdf_result.score for mdf in mdfs]) # old-mdf.magnitude
	#feasibility_point = np.interp(0,mdfs,conc_range)
	#print('MDF = 0 at approximately ' + str(feasibility_point) + ' M')
	#maximum_MDF = np.amax(mdfs)
	#print('Maximum MDF = ' + str(maximum_MDF) + 'kJ/mol')
	#flattening_point = np.interp(maximum_MDF,mdfs,conc_range)
	#print(flattening_point)
	#fig = plt.plot(conc_range,mdfs)
	#plt.xscale("log")
	return mdfs

#Calculates the MDF across two ranges of concentrations (in Molar) of two metabolites in the pathway, returning a 2d array as a result
def MDF_double_conc_sweep(cmpd_x_acronym, cmpd_y_acronym, conc_range_x, conc_range_y, pathway_file_name, compound_settings, physiological_context, custom_dGs = True, y=True):
	cs = compound_settings
	pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x=y)
	cmpd_x = cs.acr2cmpd_dict.get(cmpd_x_acronym)
	cmpd_y = cs.acr2cmpd_dict.get(cmpd_y_acronym)
	if cmpd_x is None:
		print(cmpd_x_acronym + " not a valid acronym!")
	if cmpd_y is None:
		print(cmpd_y_acronym + " not a valid acronym!")
	x = len(conc_range_x)
	y = len(conc_range_y)
	results = np.empty((x,y), dtype = PathwayMDFData)
	mdfs = np.empty((x,y))
	nn = x * y

	for i in range(x):
		for j in range(y):
			pathway._bounds.set_bounds(cmpd_x, Q_(conc_range_x[i], "M"), Q_(conc_range_x[i], "M")) # old-SetBounds
			pathway._bounds.set_bounds(cmpd_y, Q_(conc_range_y[j], "M"), Q_(conc_range_y[j], "M"))
			mdf_result = pathway.mdf_analysis()    # old-calc_mdf()
			results[i,j] = mdf_result
			mdfs[i,j] = mdf_result.mdf.magnitude
			print( 'progress {:3.1f}'.format( (i*y+j)/(x*y)*100 ) )

	return mdfs, results

#Calculates the MDF across a range of changing ratios between two metabolites, while these maintain a same total concentration
#We distinguish between the numerator metabolite and the denominator metabolite
def MDF_ratio_sweep(num_compound_acronym, den_compound_acronym, ratio_range, total_conc,pathway_file_name, compound_settings, physiological_context, custom_dGs = True, y = True):
	if type(total_conc) is float: total_conc = Q_(total_conc,"M") #if it is given as a float we assume the concentration to be in molar
	cs = compound_settings

	if y == True:
		pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x = True)
	else:
		pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x = y)
	
	mdf_results = []
	num_compound = cs.acr2cmpd_dict.get(num_compound_acronym)
	if num_compound is None:
		print(num_compound_acronym + " not a valid acronym!")
	den_compound = cs.acr2cmpd_dict.get(den_compound_acronym)
	if den_compound is None:
		print(den_compound_acronym + " not a valid acronym!")

	den_conc_range = np.divide(total_conc,(1+ratio_range))
	num_conc_range = np.multiply(ratio_range,den_conc_range)

	for i in range(ratio_range.size):
		pathway._bounds.set_bounds(num_compound, num_conc_range[i], num_conc_range[i]) # old-SetBounds
		pathway._bounds.set_bounds(den_compound, den_conc_range[i], den_conc_range[i])
		mdf_result = pathway.mdf_analysis()     # old-calc_mdf()
		mdf_results.append(mdf_result)
	mdfs_mag = np.array([mdf_result.score for mdf_result in mdf_results]) # old-np.array([mdf.magnitude for mdf in mdfs])
	#print(mdfs_mag)
	#feasibility_point = np.interp(0,mdfs_mag,ratio_range)
	#print('MDF = 0 at approximately ' + str(feasibility_point))
	#maximum_MDF = np.amax(mdfs_mag)
	#my-(optional) Auswahlverfahren für ratio
	x = ratio_range
	val = []
	delta1 = -1000
	k = -1000
	for i,j in enumerate(mdfs_mag):
    		delta2 = delta1
    		delta1 = j-k
    		k = j
    		if delta1 > -0.0001 and delta1 < 0.0001 and delta1 != -1000 and delta2 != -1000 and k != -1000:
        		val.append(x[i-1])
    		elif delta1 > 0 and delta2 < 0 and delta1 != -1000 and delta2 != -1000 and k != -1000:
        		val.append(x[i-1])
    		elif delta1 < 0 and delta2 > 0 and delta1 != -1000 and delta2 != -1000 and k != -1000:
        		val.append(x[i-1])
	if val != []:
		print('max_MDF_ratio = [',val[0],val[-1],'] -> choose')

	#print('Maximum MDF = ' + str(maximum_MDF) + 'kJ/mol')
	#flattening_point = np.interp(maximum_MDF,mdfs_mag,ratio_range)  # old-np.interp(maximum_MDF,mdf,conc_range)
	#print('flattening_point =',flattening_point) 
	#fig = plt.plot(conc_range,mdfs)
	#plt.xscale("log")
	return mdf_results,num_conc_range,den_conc_range

#Calculates the MDF across two ranges of metabolite ratios of two pairs of metabolites in the pathway, returning a 2d array as a result
def MDF_double_ratio_sweep(ratioX_tuple, ratioY_tuple, pathway_file_name, compound_settings, physiological_context, custom_dGs = True, y = True):
	(n_acr_X, d_acr_X, rr_X, tc_X) =  ratioX_tuple
	(n_acr_Y, d_acr_Y, rr_Y, tc_Y) =  ratioY_tuple
	if type(tc_X) is float: tc_X = Q_(tc_X,"M")
	if type(tc_Y) is float: tc_Y = Q_(tc_Y,"M")
	cs = compound_settings

	if y == True:
		pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x = True)
	else:
		pathway = create_pathway(pathway_file_name,cs,physiological_context,custom_dGs=custom_dGs, x = y)
	
	x = len(rr_X)
	y = len(rr_Y)
	results = np.empty((x,y), dtype = PathwayMdfSolution)    # old- PathwayMDFData
	mdfs = np.empty((x,y))
	acronyms = [n_acr_X, d_acr_X, n_acr_Y, d_acr_Y]
	compounds = [cs.acr2cmpd_dict.get(acronym) for acronym in acronyms]
	n_cmpd_X, d_cmpd_X, n_cmpd_Y, d_cmpd_Y = compounds
	print(" not a valid acronym! ".join([acronym for (acronym, compound) in zip(acronyms, compounds) if compound is None]))
	d_cr_X = np.divide(tc_X,(1+rr_X))
	n_cr_X = np.multiply(rr_X,d_cr_X)
	d_cr_Y = np.divide(tc_Y,(1+rr_Y))
	n_cr_Y = np.multiply(rr_Y,d_cr_Y)

	for i in range(x):
		for j in range(y):
			pathway._bounds.set_bounds(n_cmpd_X, n_cr_X[i], n_cr_X[i])  # old-SetBounds
			pathway._bounds.set_bounds(d_cmpd_X, d_cr_X[i], d_cr_X[i])
			pathway._bounds.set_bounds(n_cmpd_Y, n_cr_Y[j], n_cr_Y[j])
			pathway._bounds.set_bounds(d_cmpd_Y, d_cr_Y[j], d_cr_Y[j])
			mdf_result = pathway.mdf_analysis()
			results[i,j] = mdf_result
			mdfs[i,j] = mdf_result.score	
	return mdfs, results

#Returns the indexes of metabolites with MDF shadow prices above 1e-5 from the list of compounds in the compound dictionary of the MDF results
def significant_rxn_shadowprice(mdf_result):
	sps = mdf_result.reaction_df["shadow_price"].values #all values are positive
	#indexes = [idx for idx, val in enumerate(sps) if val == max(sps)]
	indexes = [idx for idx, val in enumerate(sps) if val > 1e-5]
	return indexes

#Returns a list of the ids of the reactions with significant shadow price values
def pathway_bottlenecks(mdf_result):
	indexes = significant_rxn_shadowprice(mdf_result)
	return mdf_result.reaction_df['reaction_id'].values[indexes]

#Returns a pathway with specified reaction(s) removed
def remove_reaction(pathway, rids):
	if type(rids) is not list:
		rids = [ rids ]
	rxns = pathway.reactions
	rid2idx = dict([(rxn.rid,idx) for idx, rxn in enumerate(rxns)])
	idxs = []
	for rid in rids:
		i = rid2idx.get(rid)
		if i is None:
			print(rid + " reaction not found!")
		else:
			idxs.append(i)
	new_rxns = rxns.copy()
	flxs = pathway.fluxes.copy()
	dgs = pathway.standard_dg_primes.copy()
	for i in idxs:	
		new_rxns.pop(i)
		flxs = np.delete(flxs.magnitude, i) * flxs.units
		if dgs is not None:
			dgs = np.delete(dgs.magnitude, i) * dgs.units	
	return equilibrator_pathway.Pathway(new_rxns, flxs,standard_dg_primes=dgs,bounds=pathway._bounds,config_dict = pathway.config_dict)

#Takes all reactions with signifficant shadow prices in MDF, and calculates the resulting MDF values if each of these reaction were to be removed
def reaction_removal_analysis(pathway_file_name, cmpd_settings, physiological_context, custom_dGs = True):
	(reactions,fluxes,std_dGs) = custom_parse_reactions(pathway_file_name,cmpd_settings,custom_dGs = custom_dGs)
	
	pathway = equilibrator_pathway.Pathway(reactions, fluxes,standard_dg_primes=std_dGs,bounds=cmpd_settings.bounds,config_dict = physiological_context)

	highest_mdf = Q_(-100,"kJ/mol")
	btnck_rxn = "none"
	remove_rxn_ptwy = None
	mdf_result = pathway.calc_mdf()
	indexes = significant_rxn_shadowprice(mdf_result)
	for i in indexes:
		rxns = reactions.copy()
		rxns.pop(i)
		flxs = np.delete(fluxes.magnitude, i) * fluxes.units
		dgs = None
		if std_dGs is not None:
			dgs = np.delete(std_dGs.magnitude, i) * std_dGs.units
		ptwy = equilibrator_pathway.Pathway(rxns, flxs,standard_dg_primes=dgs,bounds=cmpd_settings.bounds,config_dict = physiological_context)
		mdf = ptwy.calc_mdf().mdf
		print("Remove " + reactions[i].rid + " => " + str(mdf))
		if mdf > highest_mdf:
			highest_mdf = mdf
			btnck_rxn = reactions[i].rid
			remove_rxn_ptwy = ptwy
	return (btnck_rxn, highest_mdf, remove_rxn_ptwy)


#Calculates efficacy of reaction (note that MDFs need to be made negative)
def efficacy(dGr,temperature):
	R = Q_(8.31e-3, "kJ / mol / K")
	x = math.exp(-dGr/(R*temperature))
	return (x-1)/(x+1)

#Calculates the reversibility factor, eta_rev, based on a specified drG and temperature
def reversibility_factor(dGr,temperature):
	R = Q_(8.31e-3, "kJ / mol / K")
	x = math.exp(dGr/(R*temperature))
	return (1-x)

#Calculates the drG necessary to obtain a specificed reversibility factor, eta_rev, at a specified temperature
def EtaRev_2_dG(eta_rev,temperature):
	R = Q_(8.31e-3, "kJ / mol / K")
	return R*temperature* math.log(1-eta_rev)

#Performs Concentration Variability Analysis of the metabolites involved in a pathway, with the minimum driving force defined, in float as kJ/mol or as a quantity.
#Specific compounds can be defined as a list of compound acronyms. When no compounds are defined, the analysis is performed for all compounds in the pathway
#Results are plotted unless print_results is set as False
def MDF_CVA(MDF, pathway_file_name,sheetname, cmpd_settings, physiological_context, compounds = None, custom_dGs = True, print_results = True, y = 0):
	if type(MDF) is Q_:
		MDF = MDF.m_as('kJ/mol')
	MDF_overRT = MDF/RT.magnitude
	pathway = create_pathway(pathway_file_name,sheetname,cmpd_settings,physiological_context,custom_dGs=custom_dGs, x=y)
	if compounds is None:
		compounds = pathway.compound_ids.copy()    #old-pathway.compound_names.copy()
	compounds.remove('H2O')
	results = []
	for cmpd in compounds:
		results.append(concentration_variability_analysis(MDF_overRT, cmpd, pathway))
	if print_results:
		fig = plt.figure(figsize=(8, 5))
		ax = plt.subplot()
		CVA_plot(results, ax=ax)
		plt.ylabel('Allowed metabolite concentration (M)',fontsize=12)
		plt.xlabel('Metabolites',fontsize=12)
		ax.set_title('Concentration Variability Analysis of ' + str(sheetname[3:]) +' for MDF = ' + str(MDF) + ' kJ/mol', fontsize=12)
		plt.tight_layout()
	return results

#Plotting function with CVA result tuple list as input
def CVA_plot(cva_results,ax = None):
	if ax is None:
		ax = plt.subplot()
	ax.set_yscale("log")
	width = 0.8    
	ind = range(len(cva_results))
	cmpds, minimum, maximum = zip(*cva_results)
	maximum = [1e6 if m == 'unbounded' else m for m in list(maximum)]
	minimum = [1e-10 if m == 'unbounded' else m for m in list(minimum)]
	height = [t-b for b, t in zip(minimum,maximum)]
	ax.bar(ind, height, width, bottom=minimum)   
	ax.set_xticks(ind)
	xticks = list(cmpds)
	ax.set_xticklabels(xticks, size="medium", rotation=90)
	ax.set_ylim(1e-6, 1e-1)

#Internal function for the linear programming optimization for a single metabolite's CVA
def concentration_variability_analysis(MDF_overRT, cmpd_acr, pathway):
	
	def get_primal_array(l):
		return np.array([v.primal for v in l], ndmin=1)

	compound_idx = pathway.compound_ids.index(cmpd_acr)    #old-pathway.compound_names.index(cmpd_acr)
	cmpd_names_except = pathway.compound_ids.copy()
	cmpd_names_except.remove(cmpd_acr)
	reaction_names = [rxn.rid for rxn in pathway.reactions]

	thermo_model = pathway
	print(thermo_model.I_dir)
	inds = np.nonzero(np.diag(thermo_model.I_dir))[0].tolist()
	#x is the vector of (y | log-conc | B)
	
	A12 = thermo_model.I_dir[inds] @ thermo_model.S.T
	A13 = np.ones((len(inds), 1))
	# log conc ub and lb - remove constraints for tested metabolite
	A32 = np.eye(thermo_model.Nc)
	A32 = np.delete(A32, compound_idx, 0)
	A33 = np.zeros((thermo_model.Nc-1, 1))

	# upper bound values
	print(thermo_model.standard_dg_primes.magnitude)
	b1 = -thermo_model.standard_dg_primes.magnitude / RT.magnitude
	
	A = np.vstack( ( np.hstack( (A12, A13) ),   # driving force
					 np.hstack( (A32, A33) ),  # log conc ub
					 np.hstack( (-A32, A33) ) ) ) # log conc lb
	
	b = np.hstack(
		[b1, np.delete(thermo_model.ln_conc_ub, compound_idx, 0), -np.delete(thermo_model.ln_conc_lb, compound_idx, 0)]
	)

	#y = [Variable("y_" + r_name) for r_name in reaction_names]

	# ln-concentration variables
	lnC = [Variable("lnC_" + c_name) for c_name in pathway.compound_ids]   #-old-last .compound_names
	B = Variable("mdf")
	x = lnC + [B]

	lp_max = Model(name="CVA_MAX")
	lp_min = Model(name="CVA_MIN")

	cnstr_names = (
		["driving_force_" + r_name for r_name in reaction_names]
		+ ["log_conc_ub_" + c_name for c_name in cmpd_names_except]
		+ ["log_conc_lb_" + c_name for c_name in cmpd_names_except]
	)
	
	constraints = []
	for j in range(A.shape[0]):
		row = [A[j, i] * x[i] for i in range(A.shape[1])]
		constraints.append(
			Constraint(sum(row), ub=b[j], name=cnstr_names[j])
		)
	
	constraints.append(Constraint(-B,ub =-MDF_overRT, name='MDF_lb'))
	
	for const in constraints:
		print(const)
	
	results = []
	
	for lp, direction in [(lp_min,'min'),(lp_max,'max')]:
		lp.add(constraints)
		lp.objective = Objective(lnC[compound_idx], direction=direction)
		
		print(lp.objective)
		print("Objective value:", lp.objective.value)
		print("")
		for var in lp.variables:
			print(var.name, ":", var.primal)
		
		status = lp.optimize()
		print( status )
		if status != "optimal": 
			print(cmpd_acr + ' ' + direction + ': ' + status)
			results.append(status) 
		else:
			results.append(np.exp(lp.objective.value))
			#print( cmpd_acr, direction, lp.objective.value )

	return (cmpd_acr,results[0],results[1])




"""
FUNCTIONS FOR ECM
"""

#Returns total enzyme volume of an ECM solution
def total_enzyme_volume(model, ln_conc):
	enz_mws = list(model.ecf.mw_enz)
	enz_conc = list( np.exp( model.ecf.ECF(ln_conc).value) )
	enz_vols = [mw*conc for (mw,conc) in zip(enz_mws,enz_conc)]
	total_vol = sum(enz_vols)
	#ordered_vols = sorted(zip(enz_vols, model.reaction_ids), reverse=True)
	#ordered_percs = [(rid, vol/total_vol) for (vol,rid) in ordered_vols]
	#print(ordered_percs)
	return sum(Q_(enz_vols,'g/L'))

#Returns dictionaries of individual enzyme demand in mass concentration, in both total values and percentage of total enzyme volume from a ECM solution
#Dictionaries are ordered from highest mass concentration to lowest
def enzyme_cost_distribution(model, ln_conc):
	enz_mws = list(model.ecf.mw_enz)
	enz_conc = list( np.exp( model.ecf.ECF(ln_conc).value) )
	enz_vols = [mw*conc for (mw,conc) in zip(enz_mws,enz_conc)]
	total_vol = sum(enz_vols)
	abs_enz_vols = (sorted(zip(model.reaction_ids,list(enz_vols)), reverse=True, key=lambda x: x[1]))
	rel_enz_vols = [(rid,vol/total_vol) for (rid,vol) in abs_enz_vols]
	abs_enz_vols = dict(abs_enz_vols)
	rel_enz_vols = dict(rel_enz_vols)
	return (abs_enz_vols, rel_enz_vols)

#Returns dictionaries of enzyme demand in mass concentration distributed in cost partitions and enzyme, in both total values and percentage of total enzyme volume from a ECM solution
#Dictionaries are ordered from highest mass concentration to lowest
#ONLY FOR ECF3!!!!
def enzyme_partition_distribution(model, ln_conc): 
	partitions = model.ecf.get_enzyme_cost_partitions(ln_conc)
	distrib = []
	names = ['Capacity', 'Thermodynamic', 'Saturation', 'Sinergy']
	#for row, rid, mw in (partitions, model.reaction_ids, model.ecf.mw_enz):
	for i in range(len(model.reaction_ids)):
		ecf1 = partitions[i][0]
		rev_fraction = partitions[i][1] - 1
		kin_fraction = partitions[i][2] - 1
		enz_concs = [ecf1, ecf1*rev_fraction, ecf1*kin_fraction, ecf1*rev_fraction*kin_fraction]
		distrib.extend([(conc*model.ecf.mw_enz[i], model.reaction_ids[i], name) for conc, name in zip(enz_concs, names)])
		print()
	abs_enz_vols = sorted(distrib, reverse=True, key=lambda x: x[0])
	total_vol = total_enzyme_volume(model,ln_conc).magnitude
	rel_enz_vols = [(vol/total_vol, rid, name) for (vol, rid, name) in abs_enz_vols]
	return (abs_enz_vols, rel_enz_vols)

#Function to check if values added to tuple list are actually values
def isfloat(value):
	try:
		float(value)
		return True
	except ValueError:
		return False

#Extracting multiple kM from a single cell in csv file
def kM_dictionary(reaction_list, kM_list):
	kM_sets = []
	for kM_set_of_values in kM_list:
		#print(kM_set_of_values)
		kM_sets.append(dict([tuple(kM_pair.split(":")) for kM_pair in kM_set_of_values.split()]))
	return dict(zip(reaction_list,kM_sets))

#Check whether kinetic parameters of enzymes respect the Haldane Relationship (printing the fold difference of those which don't)
#Tolerance can be set to different % values of max difference
def check_parameters(pathway_file_name, cmpd_settings, physiological_context,abs_flux,tolerance = 0,custom_dGs=True):
	model = ECM_analysis(pathway_file_name, cmpd_settings, physiological_context,abs_flux,custom_dGs=custom_dGs)
	
	with open(pathway_file_name + '.csv','r') as csv_file:
		enzyme_data = pd.read_csv(csv_file)
		csv_file.close()

	kcatf_log = np.array([np.log(float(kcat)) for kcat in enzyme_data["kcrf(1/s)"]])
	kcatr_log = np.array([np.log(float(kcat)) for kcat in enzyme_data["kcrr(1/s)"]])

	D_S_coeff = np.diag( model.ecf.S_subs.T @ np.log(model.ecf.KMM))
	D_P_coeff = np.diag( model.ecf.S_prod.T @ np.log(model.ecf.KMM))
	
	rhs = kcatf_log + D_P_coeff - kcatr_log - D_S_coeff
	lhs = -model.ecf.standard_dg_over_rt
	ln_dif = rhs - lhs
	exp_dif = np.exp(ln_dif)
	Keq = np.exp(lhs)
	for rid, dif in zip(model.reaction_ids, exp_dif):
		if dif > (1+tolerance) or dif < (1-tolerance):
			print('In ' + rid + ', parameters yield ' + str(dif) + ' times the Keq')
	return Keq, exp_dif

#Function to prepare all data necessary for ECM analysis into required data frame
def obtain_enzyme_parameters(filename,pathway,compound_settings):
	cols = ["QuantityType", "Value", "Compound", "Reaction", "Unit"]
	list_tuples = []
	
	#Getting kEQ (equilbrium constant)
	#R = Q_(8.31e-3, "kJ / mol / K")
	#dGrs = list(pathway.standard_dg_primes)
	#kEQs = [str(math.exp(-dGr/(R*pathway.temperature))) for dGr in dGrs]
	#rks = zip(kEQs,list(pathway.reaction_ids))
	#list_tuples.extend([("equilibrium constant", kEQ, "", rid,"dimensionless") for (kEQ, rid) in rks])
	#print([("equilibrium constant", kEQ, "", rid,"dimensionless") for (kEQ, rid) in rks])
	
	with open(filename + '.csv','r') as csv_file:
		enzyme_data = pd.read_csv(csv_file)
		csv_file.close()
		
	#Getting kC (catalytic rate constant geometric mean)
	#rks = zip(enzyme_data["kC(1/s)"] , enzyme_data["Reaction Name"])
	#list_tuples.extend([("catalytic rate constant geometric mean", kC, "", rid,"1/s") for (kC, rid) in rks])
			  
	#Getting kM (Michaelis constant) -> compound and reaction specific!
	kM_dict = kM_dictionary(enzyme_data["Reaction Name"],enzyme_data["kM(mM)"])
	for rxn in pathway.reactions:
		for compound in rxn.keys():
			acronym = compound_settings.cmpd2acr_dict[compound]
			kM = kM_dict[rxn.rid].get(acronym)
			if (kM is None):
				if(acronym != "H2O"):
					print("kM for " + acronym + " was not found in " + rxn.rid)
			else:
				list_tuples.append(("Michaelis constant",kM,acronym,rxn.rid,"mM"))
	
	#Getting kcrf (catalytic rate constant of forward reaction)
	rks = zip(enzyme_data["kcrf(1/s)"] , enzyme_data["Reaction Name"])
	list_tuples.extend([("substrate catalytic rate constant", kcrf, "", rid,"1/s") for (kcrf, rid) in rks])
			  
	#Getting kcrr (catalytic rate constant of reverse reaction)
	rks = zip(enzyme_data["kcrr(1/s)"] , enzyme_data["Reaction Name"])
	list_tuples.extend([("product catalytic rate constant", kcrr, "", rid,"1/s") for (kcrr, rid) in rks])
			  
	#Getting MWe(Protein Molecular mass of enzymes)
	rks = zip(enzyme_data["MWe(Da)"] , enzyme_data["Reaction Name"])
	list_tuples.extend([("protein molecular mass", MWe, "", rid,"Da") for (MWe, rid) in rks])
			  
	#Getting MWc(Compound Molecular mass)
	compounds = Pathway.get_compounds(pathway.reactions)
	#print([("molecular mass", str(compound.mass), compound_settings.cmpd2acr_dict[compound],"","Da")])
	list_tuples.extend([("molecular mass", str(compound.mass), compound_settings.cmpd2acr_dict[compound],"","Da") for compound in compounds])

	#list_tuples = [(i,j,k,l,m) for (i,j,k,l,m) in list_tuples if isfloat(j)]
	#print([(i,j,k,l,m) for (i,j,k,l,m) in list_tuples if j == None])
	return pd.DataFrame(list_tuples, columns = cols)

#Returns model object for ECM analysis based on pathway and conditions
def ECM_analysis(pathway_file_name, cmpd_settings, physiological_context,abs_flux,custom_dGs=True):
	cs = cmpd_settings
	#print(cs.bounds.get_bound_tuple(cs.acr2cmpd_dict['2HB']))         # old-GetBoundTuple
	pathway = create_pathway(pathway_file_name,cs,physiological_context,absolute_flux=abs_flux,custom_dGs=custom_dGs)
	parameter_df = obtain_enzyme_parameters(pathway_file_name,pathway,cs)
	#print(pathway.fluxes)
	return EnzymeCostModel(pathway,parameter_df)

#Performs a Michaellis Constant parameter sweep of a specified enzyme and one of its substrate/product
def ECM_kM_sweep(kM_tuple,pathway_file_name,cs,pc,abs_flux,n_iter = 100,custom_dGs=True):
	pathway = create_pathway(pathway_file_name,cs,pc,absolute_flux=abs_flux,custom_dGs=custom_dGs)

	p_df = obtain_enzyme_parameters(pathway_file_name,pathway,cs)
	(compound_acronym,reaction_acronym,value_range) = kM_tuple
	r_idx = np.where((p_df.QuantityType == 'Michaelis constant') & (p_df.Reaction == reaction_acronym) & (p_df.Compound == compound_acronym))
	r_idx = r_idx[0][0]
	#cmpd_idx = new_model.compound_ids.index(compound_acronym)
	#rxn_idx = new_model.reaction_ids.index(reaction_acronym)
	#kM = Q_(float(p_df.loc[r_idx,'Value']),p_df.loc[r_idx,'Unit'])
	kM = float(p_df.loc[r_idx,'Value'])
	print('kM = ' + str(kM) + ' ' + p_df.loc[r_idx,'Unit'])
	parameter_range = np.multiply(kM,value_range)
	tevs = []
	lnCs = []
	new_model = EnzymeCostModel(pathway,p_df)

	for value in parameter_range:
		p_df.loc[r_idx,'Value'] = value
		new_model = EnzymeCostModel(pathway,p_df)
		lnC = new_model.ecf.optimize_ecm()[1]
		tevs.append(total_enzyme_volume(new_model, lnC).magnitude)
		lnCs.append(lnC)
	return np.array(tevs), lnCs

#This function is almost a carbon-copy of PlotEnzymeDemandBreakdown in ecm_model.py, it was edited because the original
#script would yield an error if no measured concentrations were given.
def ECM_plot(model,ln_conc,ax,fontsize_ylabel = 12):
	top_level = 3
	
	assert top_level in range(1, 5)
	costs = model.ecf.get_enzyme_cost_partitions(ln_conc)
	#print(costs)
	# give all reactions with zero cost a base value, which we will
	# also set as the bottom ylim, which will simulate a "minus infinity"
	# when we plot it in log-scale
	base = min(filter(None, costs[:, 0])) / 2.0
	idx_zero = costs[:, 0] == 0
	costs[idx_zero, 0] = base
	costs[idx_zero, 1:] = 1.0
	bottoms = np.hstack(
		[np.ones((costs.shape[0], 1)) * base, np.cumprod(costs, 1)]
		)
	steps = np.diff(bottoms)

	labels = EnzymeCostFunction.ECF_LEVEL_NAMES[0:top_level]

	ind = range(costs.shape[0])  # the x locations for the groups
	width = 0.8
	ax.set_yscale("log")
	colors = ["tab:blue", "tab:orange", "tab:brown"]
	for i, label in enumerate(labels):
		ax.bar(
			ind,
			steps[:, i].flat,
			width,
			bottom=bottoms[:, i].flat,
			color=colors[i],
			label=label,
		)

	ax.set_xticks(ind)
	xticks = model.reaction_ids
	ax.set_xticklabels(xticks, size="medium", rotation=90)
	ax.legend(loc="best", framealpha=0.2)
	ax.set_ylabel("enzyme demand [M]",fontsize = fontsize_ylabel)
	ax.set_ylim(bottom=base)

#This function is an attempt to recreate the costfunction's optimization algorithm to allow for better control
#of the optimization
def myECM(ecf: EnzymeCostFunction, method: str, ln_conc0: np.ndarray = None, n_iter: int = 10) -> np.ndarray:
		METABOLITE_WEIGHT_CORRECTION_FACTOR = 1
		"""
			Use convex optimization to find the y with the minimal total
			enzyme cost per flux, i.e. sum(ECF(lnC))
		"""

		def optfun(ln_conc: np.ndarray) -> float:
			"""
				regularization function:
					d      = x - 0.5 * (x_min + x_max)
					lambda = median(enzyme cost weights)
					reg    = 0.01 * lambda * 0.5 * (d.T * d)
			"""
			ln_conc = np.array(ln_conc.flatten(), ndmin=2).T

			# if some reaction is not feasible, give a large penalty
			# proportional to the negative driving force.
			minimal_df = ecf._DrivingForce(ln_conc).min()
			if minimal_df <= 0:
				return 1e20 * abs(minimal_df)

			enz_conc = np.exp( ecf.ECF(ln_conc) )
			met_conc = np.exp(ln_conc)

			e = float(enz_conc.T @ ecf.mw_enz)
			m = float(met_conc.T @ ecf.mw_met)
			if np.isnan(e) or e <= 0:
				raise Exception(
					"ECF returns NaN although all reactions are feasible"
				)

			if (
				ecf.regularization is None
				or ecf.regularization.lower() == "none"
			):
				return e
			elif ecf.regularization.lower() == "volume":
				return e + METABOLITE_WEIGHT_CORRECTION_FACTOR * m
			elif ecf.regularization.lower() == "quadratic":
				d = ln_conc - 0.5 * (ln_conc.min() + ln_conc.max())
				return e + QUAD_REGULARIZATION_COEFF * 0.5 * float(d.T * d)
			else:
				raise Exception(
					"Unknown regularization: " + ecf.regularization
				)

		if ln_conc0 is None:
			mdf, ln_conc0 = ecf.MDF()
			if np.isnan(mdf) or mdf < 0.0:
				raise ValueError(
					"It seems that the problem is thermodynamically"
					" infeasible, therefore ECM is not applicable."
				)
		assert ln_conc0.shape == (ecf.Nc, 1)

		bounds = list(zip(ecf.ln_conc_lb.flat, ecf.ln_conc_ub.flat))

		min_res = np.inf
		ln_conc_min = ln_conc0
		for i in range(n_iter):
			ln_conc0_rand = ln_conc0 * (
				1.0 + 0.1 * np.random.rand(ln_conc0.shape[0], 1)
			)
			r = minimize(
				optfun, x0=ln_conc0_rand, bounds=bounds, method=method#'SLSQP'#'L-BFGS-B'
			)

			if not r.success:
				logging.info(
					f"iteration #{i}: optimization unsuccessful "
					f"because of {r.message}, "
					"trying again"
				)
				continue

			res = optfun(r.x)
			if res < min_res:
				if min_res == np.inf:
					logging.info("iteration #%d: cost = %.5f" % (i, res))
				else:
					logging.info(
						"iteration #%d: cost = %.5f, decrease factor = %.3e"
						% (i, res, 1.0 - res / min_res)
					)
				min_res = res
				ln_conc_min = np.array(r.x, ndmin=2).T
			else:
				logging.info(
					"iteration #%d: cost = %.5f, no improvement" % (i, res)
				)

		return ln_conc_min

#Misc functions

#Return index of first element above a certain value from an ordered 1d array
def idx_of_val_just_above(sorted_array, target_value):
	i= 0
	if target_value > sorted_array[-1]:
		print('Value ' + str(target_value) + ' too big for array!')
		return -1
	while sorted_array[i] < target_value:
		i += 1
	return i
	
