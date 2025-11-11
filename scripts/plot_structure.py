from Bio.PDB import MMCIFParser, DSSP
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import PathPatch
import matplotlib.path as mpath
import matplotlib.patches as patches
import matplotlib.pylab as pylab
import re

pylab.rcParams['pdf.fonttype'] = 42
pylab.rcParams['ps.fonttype'] = 42

def merge_intervals(nums):
    if not nums:
        return []

    nums = sorted(nums)
    intervals = []
    start = nums[0]
    end = nums[0]

    for i in range(1, len(nums)):
        if nums[i] == end + 1:
            end = nums[i]
        else:
            intervals.append((start, end))
            start = nums[i]
            end = nums[i]

    intervals.append((start, end))
    return intervals

def get_structure(pbd_file,
                  map: dict = None):

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", pbd_file) 

    model = structure[0] 
    dssp = DSSP(model, pbd_file)

    mmcif_dict = MMCIF2Dict(pbd_file)
    entity_names = mmcif_dict["_entity.pdbx_description"]
    chain_ids = mmcif_dict["_struct_asym.id"]
    entity_chain_map = mmcif_dict["_struct_asym.entity_id"]

    chain_entity_map = {chain_id: entity_names[int(entity_id) - 1] for chain_id, entity_id in zip(chain_ids, entity_chain_map)}
    # print(chain_entity_map)

    histone_structure = {}

    for key in dssp.keys():
        chain_id = key[0]
        if map:
            chain_id = map[chain_id]

        residue_position = key[1][1]
        secondary_structure = dssp[key][2]
        
        entity_name = chain_entity_map.get(chain_id, "Unknown")
        histone_family = re.search(r'(H1|H2A|H2B|H3|H4)', entity_name)

        if not histone_family:
            continue

        histone_name = histone_family.group(1) + f"_{chain_id}"

        if histone_name not in histone_structure:
            histone_structure[histone_name] = {
                "helix": [],
                "sheet": []
            }

        if secondary_structure == "H":
            histone_structure[histone_name]["helix"].append(residue_position)
        elif secondary_structure == "E":
            histone_structure[histone_name]["sheet"].append(residue_position)

    for histone_name, structure_info in histone_structure.items():
        for key in structure_info:
            structure_info[key] = merge_intervals(structure_info[key])
    
    histone_structure = dict(sorted(histone_structure.items(), key=lambda x: x[0]))
    return histone_structure


def process_sequence(seq, acidic, basic):
    """Identify properties of amino acids in the sequence."""
    aa_properties = []
    for aa in seq:
        if aa in acidic:
            aa_properties.append('acidic')
        elif aa in basic:
            aa_properties.append('basic')
        else:
            aa_properties.append('neutral')
    return aa_properties

def plot_sequence(ax_seq, seq, aa_properties, per_amino_acid_width, property_colors, seq_label, neutral_color):
    """Plot the main amino acid sequence."""
    seq_len = len(seq)
    for i, (aa, prop) in enumerate(zip(seq, aa_properties)):
        x_pos = i * per_amino_acid_width
        color = property_colors.get(prop, neutral_color)
        rect = plt.Rectangle((x_pos, 0), per_amino_acid_width, 1,
                             facecolor=color, edgecolor='black',
                             linewidth=0.5, linestyle='-')
        ax_seq.add_patch(rect)
        text_color = 'black' if prop == 'neutral' else 'white'
        ax_seq.text(x_pos + per_amino_acid_width / 2, 0.5, aa, ha='center', va='center',
                    color=text_color, fontfamily='Courier')

        if i % 10 == 0:
            # Add position labels
            ax_seq.text(x_pos + per_amino_acid_width / 2, 1.5, str(i+1), ha='center', va='center', fontfamily='Courier')

    # Set axis limits and labels
    ax_seq.set_xlim(0, seq_len * per_amino_acid_width)
    ax_seq.set_ylim(-0.6, 1)
    ax_seq.set_xticks([])
    ax_seq.set_yticks([])
    ax_seq.set_aspect('auto')
    for pos in ("left", "right", "top", "bottom"):
        ax_seq.spines[pos].set_visible(False)
    
    if seq_label:
        ax_seq.text(-0.1, 0.5, seq_label, ha='right', va='center')

def add_annotations(ax_seq, annotations, per_amino_acid_width):
    """Add annotations to the sequence plot."""
    y_offset_more_lst = []
    for annotation in annotations:
        label = annotation.get('label')
        start, end = annotation.get('range')
        end += 1  # Adjust for inclusive range
        style = annotation.get('style', 'rectangle')
        width = end - start
        y_offset_more = annotation.get('y_offset_more', 0)
        y_offset_more_lst.append(y_offset_more)
        style_color = annotation.get('color', 'black')
        text_offset = annotation.get('text_offset', 0)
        # Scale positions
        start_pos = start * per_amino_acid_width
        width_scaled = width * per_amino_acid_width
        text_offset_scaled = text_offset * per_amino_acid_width

        if style == 'rectangle':
            arrow_y = -0.3 - y_offset_more
            ax_seq.annotate('', 
                            xy=(start_pos + width_scaled, arrow_y),
                            xytext=(start_pos, arrow_y), 
                            arrowprops=dict(arrowstyle='|-|,widthA=0.2, angleA=0, widthB=0.2, angleB=0', color=style_color, linewidth=0.5))

        elif style == 'mark':
            ax_seq.annotate('', xy=(start_pos + width_scaled / 2, -0.1 - y_offset_more),
                            xytext=(start_pos + width_scaled / 2, -0.3 - y_offset_more),
                            arrowprops=dict(arrowstyle='-|>,head_width=0.3,head_length=0.3', color=style_color, linewidth=0.5))
            
        elif style == 'arrow':
            arrow_y = -0.3 - y_offset_more 
            ax_seq.annotate('', 
                            xy=(start_pos + width_scaled, arrow_y),
                            xytext=(start_pos, arrow_y), 
                            arrowprops=dict(arrowstyle='->', color=style_color, linewidth=0.5))

        elif style == 'brace':
            Path = mpath.Path
            middle = start_pos + width_scaled / 2
            verts = [
                (start_pos, -0.15 - y_offset_more),
                (start_pos, -0.45 - y_offset_more),
                (middle - 0.02 * width_scaled , -0.45 - y_offset_more), 
                (middle, -0.55 - y_offset_more),
                (middle + 0.02 * width_scaled, -0.45 - y_offset_more), 
                (start_pos + width_scaled, -0.45 - y_offset_more),
                (start_pos + width_scaled, -0.15 - y_offset_more),
            ]
            codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO]
            path = mpath.Path(verts, codes)
            patch = PathPatch(path, fill=False, color=style_color, linewidth=0.5)
            ax_seq.add_patch(patch)

        elif style == 'helix':
            num_turns = width * 4 
            x_values = np.linspace(start_pos, start_pos + width_scaled, int(width * 100))
            y_offset = -0.3 - y_offset_more
            amplitude = 0.1 
            y_values = amplitude * np.sin(2 * np.pi * num_turns * (x_values - start_pos) / width_scaled) + y_offset
            ax_seq.plot(x_values, y_values, color=style_color, linewidth=0.5)

        # Add annotation text below
        ax_seq.text(start_pos + width_scaled / 2 + text_offset_scaled, -0.55 - y_offset_more, label,
                    ha='center', va='top')
    max_y_offset = max(y_offset_more_lst) if y_offset_more_lst else 0
    return max_y_offset

def plot_additional_sequences(ax_seq, add_seqs, per_amino_acid_width, seq_offset, cut_start, gap_lst=None):
    """Plot additional sequences below the main sequence."""
    seq_offset_more = 0

    for add_seq_info in add_seqs:

        seq = add_seq_info.get("seq")
        colors = add_seq_info.get("color")
        if cut_start:
            seq = seq[1:]
            colors = colors[1:]

        add_seq_label = add_seq_info.get("label")
        for i, (aa, color) in enumerate(zip(seq, colors)):

            left_color = colors[i - 1] if i > 0 else None
            right_color = colors[i + 1] if i < len(colors) - 1 else None
            x_pos = i * per_amino_acid_width

            left_equal = left_color == color if left_color else False
            right_equal = right_color == color if right_color else False

            if isinstance(color, str):
                rect = plt.Rectangle((x_pos, seq_offset + seq_offset_more), per_amino_acid_width, 1,
                                      facecolor=color, edgecolor='black',
                                      linewidth=0, linestyle='-')
            elif isinstance(color, list):
                
                bottom = seq_offset + seq_offset_more
                for cl in color:
                    rect = plt.Rectangle((x_pos, bottom), per_amino_acid_width, 1 / len(color),
                                         facecolor=cl, edgecolor='black',
                                         linewidth=0.5, linestyle='-')
                    ax_seq.add_patch(rect)
                    bottom += 1 / len(color)

            if not left_equal:
                # Add left border
                ax_seq.add_line(plt.Line2D(
                    [x_pos, x_pos],
                    [seq_offset + seq_offset_more, seq_offset + seq_offset_more + 1],
                    color='black',
                    linewidth=0.5
                ))
            if not right_equal:
                # Add right border
                ax_seq.add_line(plt.Line2D(
                    [x_pos + per_amino_acid_width, x_pos + per_amino_acid_width],
                    [seq_offset + seq_offset_more, seq_offset + seq_offset_more + 1],
                    color='black',
                    linewidth=0.5
                ))

            ax_seq.add_patch(rect)
            ax_seq.text(x_pos + per_amino_acid_width / 2, seq_offset + seq_offset_more + 0.5,
                        aa, ha='center', va='center', color="black", fontfamily='Courier')
            
        # Add upper and lower borders
        ax_seq.add_line(plt.Line2D(
            [0, len(seq) * per_amino_acid_width],
            [seq_offset + seq_offset_more, seq_offset + seq_offset_more],
            color='black',
            linewidth=0.5
        ))
        ax_seq.add_line(plt.Line2D(
            [0, len(seq) * per_amino_acid_width],
            [seq_offset + seq_offset_more + 1, seq_offset + seq_offset_more + 1],
            color='black',
            linewidth=0.5
        ))

        ax_seq.text(-0.1, seq_offset + seq_offset_more + 0.5, add_seq_label, ha='right', va='center')
        
        if not gap_lst:
            seq_offset_more -= 1
        else:
            try:
                seq_offset_more -= gap_lst.pop(0)
            except IndexError:
                seq_offset_more -= 1

    return seq_offset_more

def plot_bar_graph(ax_bar, x, values, per_amino_acid_width, bar_color, group_values, group_colors, stacked, bar_width_ratio=0.8,):
    """Plot the bar graph above the sequence."""
    if group_values is not None:
        groups = list(values.keys())
        n_groups = len(groups)
        if not stacked:
            # Side-by-side bars
            total_bar_width = bar_width_ratio * per_amino_acid_width
            bar_width = total_bar_width / n_groups

            for i, group in enumerate(groups):
                y_values = values[group]
                x_positions = x + per_amino_acid_width / 2 - total_bar_width / 2 + i * bar_width + bar_width / 2
                bar_color = group_colors.get(group, None) if group_colors else None
                ax_bar.bar(x_positions, y_values,
                           width=bar_width, label=group, color=bar_color,
                           edgecolor='black', linewidth=0.5)
        else:
            # Stacked bars
            bottom_values = np.zeros(len(x))
            for group in groups:
                y_values = values[group]
                bar_color = group_colors.get(group, None) if group_colors else None
                ax_bar.bar(x + per_amino_acid_width / 2, y_values, width=bar_width_ratio * per_amino_acid_width,
                           bottom=bottom_values, label=group, color=bar_color,
                           edgecolor='black', linewidth=0.5)
                bottom_values += y_values
    else:
        # Single group
        ax_bar.bar(x + per_amino_acid_width / 2, values, width=bar_width_ratio * per_amino_acid_width,
                   color=bar_color, align='center',
                   edgecolor='black', linewidth=0.5)

def add_axvspans(ax_bar, axvspans, per_amino_acid_width, cut_start):
    """Add vertical spans to the bar plot."""
    if axvspans is not None:
        axvspan_label_set = set()
        for axvspan_item in axvspans:
            axvspan_label = axvspan_item.get("label", None)
            axvspan_sites = axvspan_item.get("sites")
            axvspan_color = axvspan_item.get("color", "gray")
            axvspan_alpha = axvspan_item.get("alpha", 0.25)

            for axvspan_site in axvspan_sites:
                if cut_start:
                    axvspan_site -= 1
                
                if axvspan_label not in axvspan_label_set:
                    ax_bar.axvspan(axvspan_site * per_amino_acid_width,
                                   (axvspan_site + 1) * per_amino_acid_width,
                                   color=axvspan_color, alpha=axvspan_alpha, label=axvspan_label, linewidth=0)
                    axvspan_label_set.add(axvspan_label)
                else:
                    ax_bar.axvspan(axvspan_site * per_amino_acid_width,
                                   (axvspan_site + 1) * per_amino_acid_width,
                                   color=axvspan_color, alpha=axvspan_alpha, linewidth=0)
        ax_bar.legend()

def add_modification_markers(ax_bar, modification_sites, modification_colors, per_amino_acid_width, cut_start):
    """Add modification markers under the bar plot without altering y-axis limits."""

    min_y, max_y = ax_bar.get_ylim()
    ylim_range = max_y - min_y

    y_pos = min_y - per_amino_acid_width / 4 - ylim_range * 0.02

    for mod_site in modification_sites:
        site = mod_site['site']
        modifications = mod_site['modification']
        
        # Adjust site index if cut_start is True
        if cut_start:
            site -= 1  # Assuming site is 0-based index, adjust accordingly if needed
        
        x_pos = site * per_amino_acid_width + per_amino_acid_width / 2
        
        num_mods = len(modifications)
        theta1 = 0
        for mod in modifications:
            theta2 = theta1 + 360 / num_mods
            color = modification_colors.get(mod, 'gray')
            wedge = patches.Wedge((x_pos, y_pos), per_amino_acid_width / 4, theta1, theta2, facecolor=color, edgecolor='black', lw=0.5)
            ax_bar.add_patch(wedge)
            theta1 = theta2
    
    ax_bar.set_ylim(y_pos - per_amino_acid_width / 2, max_y)

    ticks = ax_bar.get_yticks()
    filtered_ticks = [tick for tick in ticks if tick >= min_y]
    ax_bar.set_yticks(filtered_ticks)
    ax_bar.tick_params(axis='y', which='minor', left=False)

    ax_bar.spines['left'].set_bounds(min_y, max_y)


def plot_peptide_sequence(
        seq: str,
        seq_label: str = None,
        seq_offset: float = -2,
        annotations: list = None,
        height_ratios: list = [4, 1],
        figheight: float = 3,
        per_amino_acid_width: float = 0.12,
        global_fontsize: int = 7,
        neutral_color: str = 'lightgray',
        acidic_color: str = 'red',
        basic_color: str = 'blue',
        values=None,
        bar_color: str = 'gray',
        bar_width_ratio: float = 0.8,
        bar_label: str = None,
        stacked: bool = False,
        hspace: float = 0.2,
        cut_start: bool = False,
        add_hline: float = None,
        add_seqs: list = None,
        gap_lst: list = None,
        group_values: dict = None,
        group_colors: dict = None,
        axvspans: list = None,
        save_fig: str = None,
        # modification_sites: list = None,
        # modification_colors: dict = None,
        ):
    """Plot a peptide sequence with optional annotations and bar graphs."""

    pylab.rcParams.update({
        'font.family': 'Arial',
        'font.size': global_fontsize,   
        'axes.titlesize': global_fontsize,      
        'axes.labelsize': global_fontsize,             
        'xtick.labelsize': global_fontsize,        
        'ytick.labelsize': global_fontsize,          
        'legend.fontsize': global_fontsize,           
        'figure.titlesize': global_fontsize,        
    })

    acidic = ['D', 'E']
    basic = ['K', 'R', 'H']

    if group_values:
        values = group_values
    
    # Adjust data if cutting the first amino acid
    if cut_start:
        seq = seq[1:]
        if group_values:
            values = {group: value[1:] for group, value in values.items()}
        else:
            values = values[1:]
        if annotations:
            adjusted_annotations = []
            for annotation in annotations:
                start, end = annotation.get('range')
                new_start = max(start - 1, 0)
                new_end = end - 1
                new_annotation = annotation.copy()
                new_annotation['range'] = (new_start, new_end)
                adjusted_annotations.append(new_annotation)
            annotations = adjusted_annotations

    # Process the sequence
    aa_properties = process_sequence(seq, acidic, basic)
    property_colors = {'acidic': acidic_color, 'basic': basic_color, 'neutral': neutral_color}
    seq_len = len(seq)
    x = np.arange(seq_len) * per_amino_acid_width  # Scale x

    # Create figure and subplots
    fig, (ax_bar, ax_seq) = plt.subplots(2, 1,
                                         figsize=(seq_len * per_amino_acid_width, figheight),
                                         sharex=True,
                                         gridspec_kw={'height_ratios': height_ratios})

    # Adjust subplot spacing
    plt.subplots_adjust(hspace=hspace)

    # Plot amino acid sequence
    plot_sequence(ax_seq, seq, aa_properties, per_amino_acid_width, property_colors, seq_label, neutral_color)

    # Add annotations
    if annotations:
        max_y_offset = add_annotations(ax_seq, annotations, per_amino_acid_width)
        ax_seq.set_ylim(-0.6 - max_y_offset, 1)
    else:
        ax_seq.set_ylim(-0.6, 1)
    
    # Plot additional sequences
    if add_seqs is not None:
        seq_offset_more = plot_additional_sequences(ax_seq, add_seqs, per_amino_acid_width, seq_offset, cut_start, gap_lst)
        ax_seq.set_ylim(-0.6 + (seq_offset + seq_offset_more), 1)

    # Plot bar graph
    if values is not None:
        # Add vertical spans if provided
        add_axvspans(ax_bar, axvspans, per_amino_acid_width, cut_start)
        plot_bar_graph(ax_bar, x, values, per_amino_acid_width, bar_color, group_values, group_colors, stacked, bar_width_ratio=bar_width_ratio)

        ax_bar.set_xticks([])
        for pos in ("right", "top", "bottom"):
            ax_bar.spines[pos].set_visible(False)
        
        ax_bar.set_ylabel(bar_label or 'Values')
        
        if add_hline:
            ax_bar.axhline(add_hline, color='red', linestyle='--', linewidth=0.5)
    else:
        ax_bar.axis('off')

    # # Add modification markers
    # if modification_sites and modification_colors:
    #     add_modification_markers(ax_bar, modification_sites, modification_colors, per_amino_acid_width, cut_start)
    #     # No need to adjust y-limits anymore

    if save_fig:
        plt.savefig(save_fig, bbox_inches='tight')
    plt.show()

    return fig, ax_bar, ax_seq