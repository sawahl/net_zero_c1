import numpy as np
from equilibrator_api import Q_, ureg
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection


# converts mdf_result data from print_MDF_results_individual/object sequence to a array element
def result_to_array(x):
    val = x.values
    ary = np.zeros(x.size)
    for i in range(x.size):
        ary[i] = np.array(val[i])
    return ary


# plot bar chart for the delta-G profile
def plot_bar_chart( title,MDF,g3,g4,r_id,shp_r, ax = None):
    
    dG_phy = result_to_array(g3)
    dG_opt = result_to_array(g4)
    reaction_id = r_id.values

    shp = result_to_array(shp_r)
    bn_ind = []
    for i,j in enumerate(shp):
        if j != 0:
            bn_ind.append(i)

    reaction_id_bn = reaction_id[bn_ind]
    dG_opt_bn = dG_opt[bn_ind]

    
    if ax == None:
        plt.figure( figsize = ( 10, 5 ), dpi = 300 )
        ax = plt.subplot()
        
    ax.bar(reaction_id, dG_phy, color = 'royalblue', label='Physiological concentrations (1 mM)')
    ax.bar(reaction_id, dG_opt, width=0.5, color = 'burlywood', label='MDF-optimized concentrations')
    ax.bar(reaction_id_bn, dG_opt_bn, width=0.5, color = 'red', label='Bottleneck reactions')
    ax.set_ylabel('$\Delta_r G^\prime$ / kJ/mol',fontsize=10)
    ax.set_title(title[3:] + ' : ' + 'MDF = {:0.2f} kJ/mol'.format(MDF),fontsize=11)
    ax.tick_params( axis='x', rotation=90 )
    ax.legend(fontsize=10)
    ax.grid(linewidth = 0.5)
	


def plot_concentrations_my(mdf_result):
        """Plot compound concentrations.

        :return: matplotlib Figure
        """
        ureg.setup_matplotlib(True)

        
        fig,ax = plt.subplots(1, 1, figsize=(10, mdf_result.compound_df.shape[0] * 0.25), dpi = 250)
        
        data_df = mdf_result.compound_df.copy()
        data_df["y"] = np.arange(0, data_df.shape[0])
        data_df["shadow_sign"] = data_df.shadow_price.apply(np.sign)

        for shadow_sign, group_df in data_df.groupby("shadow_sign"):
            if shadow_sign == -1:
                color, marker, label = ("blue", "<", "shadow price < 0")
            elif shadow_sign == 1:
                color, marker, label = ("red", ">", "shadow price > 0")
            else:
                color, marker, label = ("grey", "d", "shadow price = 0")

            group_df.plot.scatter(
                x="concentration",
                y="y",
                s=40,
                c=color,
                marker=marker,
                ax=ax,
                zorder=2,
                colorbar=False,
                label=label,
            )
        ax.set_ylabel("")
        ax.set_yticks(data_df.y)

        for j, row in enumerate(data_df.itertuples(index=True)):
            ax.plot(
                [row.lower_bound.m_as("M"), row.upper_bound.m_as("M")],
                [row.y, row.y],
                color="lightgrey",
                linewidth=3,
                zorder=1,
                label="allowed range" if j == 0 else None,
            )

        ax.set_yticklabels(data_df.compound_id, fontsize=9)
        ax.set_xlabel("Concentration / M")
        ax.set_xscale("log")
        ax.set_ylim(-0.5, data_df.shape[0] + 0.5)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        ax.grid(linewidth = 0.5)
        plt.show()
        
def plot_driving_forces_my(mdf_result):
        """Plot cumulative delta-G profiles.

        :return: matplotlib Figure
        """
        ureg.setup_matplotlib(True)

        fig, ax = plt.subplots(1, 1, figsize=(3 + mdf_result.reaction_df.shape[0] * 0.3, 5), dpi = 250)

        data_df = mdf_result.reaction_df.copy()
        data_df.reindex()

        data_df[
            "cml_dgm"
        ] = mdf_result.reaction_df.physiological_dg_prime.cumsum().apply(
            lambda x: x.m_as("kJ/mol")
        )
        data_df[
            "cml_dg_opt"
        ] = mdf_result.reaction_df.optimized_dg_prime.cumsum().apply(
            lambda x: x.m_as("kJ/mol")
        )

        xticks = 0.5 + np.arange(data_df.shape[0])
        xticklabels = data_df.reaction_id.tolist()

        yvalues_phy = [0.0] + data_df.cml_dgm.tolist()
        yvalues_opt = [0.0] + data_df.cml_dg_opt.tolist()

        ax.plot(
            yvalues_phy,
            label="Physiological concentrations (1 mM)",
            color="lightgreen",
            linestyle="--",
            zorder=1,
        )
        ax.plot(
            yvalues_opt,
            label="MDF-optimized concentrations",
            color="grey",
            linewidth=2,
            zorder=2,
        )
        lines = [
            ((i, yvalues_opt[i]), (i + 1, yvalues_opt[i + 1]))
            for i in data_df[data_df.shadow_price != 0].index
        ]
        lines = LineCollection(
            lines,
            label="Bottleneck reactions",
            linewidth=3,
            color="red",
            linestyle="-",
            zorder=3,
            alpha=1,
        )
        ax.add_collection(lines)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, rotation=45, ha="center")
        ax.set_xlim(0, data_df.shape[0])

        ax.set_ylabel(r"Cumulative $\Delta_r G^\prime$ / kJ/mol", fontsize=14)
        ax.legend(loc="best")
        ax.set_title(f"MDF = {mdf_result.score:.2f} kJ/mol",fontsize=13)
        ax.grid(linewidth = 0.5)
        plt.show()