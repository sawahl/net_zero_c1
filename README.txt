### NET ZERO C1 ###

### INTRODUCTION ###

The transition to a sustainable, circular economy will require the use of renewable feedstocks and bio-based processes.
Engineered microorganisms using renewable feedstocks like sugars from agricultural waste or other renewable carbon sources 
can operate at low process temperatures and no toxic waste is generated.Currently, industrial microorganisms use mainly sugar
based feedstocks; while microbes do grow efficiently on sugars, sugars are not optimal for the production of platform chemicals
with a high degree of reduction.  For many reduced chemicals, there is a low energy-to-carbon ratio or inefficient metabolic 
pathways starting from sugars. Consequently, only low carbon yields can be obtained and there are significant CO2 emissions.

The combination of different feedstocks with complementary properties could overcome this hurdle. Methane obtained from the
anaerobic breakdown of biomass (biogas) is abundant, and methane is the C1 molecule with the highest degree of reduction. 
These characteristics render methane and its derivatives such as methanol as great candidates for combination with sugar feedstock
to enable full carbon and energy efficiency and zero CO2 emissions.


### DESCRIPTION ###

This project is based in an in-house developed stoichiometric model coded as a JSON file. The pathways included in this model can be classified in several groups:

- Core metabolic pathways: This includes the Pentose Phosphate Pathway, the Tricarboxylic Acid Cycle, oxidative respiration and biomass generation.

- Sugar metabolizing pathways: This includes the Embden-Meyerhof-Parnas glycolysis, the Entner-Doudoroff glycolysis, and the non-oxidative glycolysis (Bogorad et al., 2013)

- Methanol incorporation pathways: Here, the included pathways can be classified in different sub-groups:
	- Natural pathways: Ribulos Monophosphate cycle (RumP) and Serine cycle
	- Synthetic serine cycle derivates: Modified serine Cycle (Yu et al., 2018), Synthetic Serine-Threonine Cycle (Yishai et al., 2017), Homoserine Cycle (He et al., 2020).
	- Synthetic  pathways based on natural reactions: Dihydroxyacetone Synthase (DAS) Pathway (De Simone et al., 2020), Reductive Glycine Pathway (Bar-Even et al., 2013)
	- Synthetic pathways based in non-natural, engineered reactions: Formolase Pathway (Siegel et al., 2015), Synthetic Acetyl-CoA (SACA) Pathway (Lu et al.,2019),
	  2-hydroxyacyl-CoA lyase (Hacl) Pathway (Chou et al., 2019), Glycolaldehyde Assimilation Pathway (GAA) Pathway (Yang et al., 2019)
          and Glycolaldehyde-allose 6-phosphate Assimilation (GAPA) Pathway (Mao et al., 2021).

- Product pathways:  pathways for the synthesis of PHB, noreugenin and butyric acid methyl esther. 

The Jupyter Notebooks in this repository allow the user to read and explore this model, as well as to run Flux Balance Analysis (FBA) and Min-Max Driving Force (MDF) using Cobrapy (Ebrahim et al., 2013) and eQuilibrator (Beber et al., 2021) respectively. Simulations of all possible combinations of sugar metabolizing and methanol incorporating pathways (up to a total of 39 different combinations) can be run with one
of the aforementioned target products as objective reaction.

The only exception is the calculation of MDFs with combinations including the Hacl Pathway, as no thermodynamic data could be retrieved for the metabolite glycolyl-CoA.

The user can modify the influx of substrates (glucose, xylose and methanol), shut down the intake/ exchange of CO2 and add desired knockouts. Finally, the user can export the data generated in the simulations in the form of an Excel file.


### EXTERNAL PYTHON PACKAGES REQUIRED ###

Package 		Version
cobra 			0.25.0
equilibrator-api	0.4.5.1
matplotlib		3.5.1
numpy			1.21.5
openpyxl		3.0.10
pandas			1.4.2
scipy			1.7.3


### HOW TO RUN ###

## FBA ##
First, download all files in a single folder. Open the Notebook Net_Zero_C1_FBA. In the second cell, select the desired version of the JSON model (either with or without enzyme coupling for PHB production). Run all consequent cells without modifying anything. The target product can be set in the cell "Set objective product": comment/uncomment the lines belonging to each of the four different products. Substrate influx can be modified in the following cell.

Simulations can be run by running the cell "Run simulation and control CO2 exchange". Inside this cell, the bounds for the exchange of CO2 can be modified; also desired knockouts can be introduced. The simulations are run using parsimonius FBA (pFBA), which optimizes the objective and minimizes the total sum of fluxes necessary. By using pFBA instead of normal FBA, we can obtain the minimum sugar to methanol ratio necessary to optimize product formation without CO2 emmisions. Once the simulation finish, a Dataframe including all fluxes for all reactions and pathway combinations will be displayed.

Finally, results can be saved in Excel files. Modify the variable filename with the desired name for your Excel file and run the cells. You will obtain an Excel file: the first sheet includes the Datafram with all simulation results; consequent sheets include detailed information of all pathway combinations which fulfill two conditions: there is product formation and the carbon efficiency (computed as moles of C in product divided by moles of C in the substrates) is higher than 99% (this value can be modified by the user).

## MDF ##

First, open the Excel sheet obtained from the FBA calculations and delete all sheets of simulations including the Hacl Pathway (otherwise an error will be raised during the MDF calculations). Change the variable results_excel with the name of this modified Excel file and run the cell "Compute MDFs". After simulations are complete, you can save the results in an Excel file by redefying the variable filename and running the cell "Save results in Excel". In addition, you can obtain plots for each simulations by running the next cells. 


### REFERENCES###

-Bogorad, I. W., Lin, T. S., & Liao, J. C. (2013). Synthetic non-oxidative glycolysis enables complete carbon conservation. Nature, 502(7473), 693-697. https://doi.org/10.1038/nature12575 
-Yu, H. and J.C. Liao, A modified serine cycle in Escherichia coli coverts methanol and CO2 to two-carbon compounds. Nature Communications, 2018. 9(1).
-Yishai, O., et al., Engineered Assimilation of Exogenous and Endogenous Formate in Escherichia coli. ACS Synth Biol, 2017. 6(9): p. 1722-1731.
-He, H., et al., An optimized methanol assimilation pathway relying on promiscuous formaldehyde-condensing aldolases in E. coli. Metab Eng, 2020. 60: p. 1-13.
-De Simone, A., et al., Mixing and matching methylotrophic enzymes to design a novel methanol utilization pathway in E. coli. Metabolic Engineering, 2020. 61: p. 315-325.
-Bar-Even, A., et al., Design and analysis of metabolic pathways supporting formatotrophic growth for electricity-dependent cultivation of microbes. Biochim Biophys Acta, 2013. 1827(8-9): p. 1039-47.
-Siegel, J.B., et al., Computational protein design enables a novel one-carbon assimilation pathway. Proc Natl Acad Sci U S A, 2015. 112(12): p. 3704-9.
-Lu, X., et al., Constructing a synthetic pathway for acetyl-coenzyme A from one-carbon through enzyme design. Nature Communications, 2019. 10(1).
-Chou, A., et al., 2-Hydroxyacyl-CoA lyase catalyzes acyloin condensation for one-carbon bioconversion. Nature Chemical Biology, 2019. 15(9): p. 900-906.
-Yang, X., et al., Systematic design and in vitro validation of novel one-carbon assimilation pathways. Metab Eng, 2019. 56: p. 142-153.
-Mao, Y., et al., Non-natural aldol reactions enable the design and construction of novel one-carbon assimilation pathways in vitro. Frontiers in Microbiology, 2021: p. 1360.
-Ebrahim, A., et al., COBRApy: COnstraints-Based Reconstruction and Analysis for Python. BMC Systems Biology, 2013. 7(1): p. 74.
-Beber, M.E., et al., eQuilibrator 3.0: a database solution for thermodynamic constant estimation. Nucleic Acids Research, 2021. 50(D1): p. D603-D609.