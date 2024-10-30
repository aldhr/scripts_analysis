import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def read_xvg(file_path):
    """Lee un archivo .xvg y genera un DataFrame de pandas."""
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith(('#', '@')):
                data.append([float(x) for x in line.split()])
    return pd.DataFrame(data, columns=['Residue', 'RMSF'])

def modify_xvg(df, shift_value=0):
    """Modifica el índice para asegurar que coincida con la
    numeración de Wuhan-Hu-1"""
    df['Residue'] += shift_value
    return df

# Rutas de archivos y configuraciones
files = {
    'Wuhan-Hu1': ('wt/Wuhan-Hu1.xvg', 0),
    'Delta21J1': ('del1/delta21j1.xvg', 2),
    'Delta21J2': ('del2/delta21j2.xvg', 2),
    'Omicron21K': ('omi21k/omicron21k.xvg', 3),
    'Omicron21L': ('omi21l/omicron21l.xvg', 3),
    'Omicron22B': ('omi22b/omicron22b.xvg', 5),
    'Omicron22E': ('omi22e/omicron22e.xvg', 6)
}

highlight_regions = {  
    (474, 488): 'yellow',
    (330, 345): 'darkorange',
    (358, 374): 'hotpink',
}

# Agregar información de mutaciones
mutations = {
    'Wuhan-Hu1': [],
    'Delta21J1': ['L452R', 'T478K'],
    'Delta21J2': ['L434V', 'L452R', 'T478K'],
    'Omicron21K': [
        'G339D', 'S371L', 'S373P', 'S375F', 'K417N', 
        'N440K', 'G446S', 'S477N', 'T478K', 'E484A', 
        'Q493R', 'G496S', 'Q498R', 'N501Y', 'Y505H'
    ],
    'Omicron21L': [
        'G339D', 'S371F', 'S373P', 'S375F', 'T376A', 
        'D405N', 'R408S', 'K417N', 'N440K', 'S477N', 
        'T478K', 'E484A', 'Q493R', 'Q498R', 'N501Y', 'Y505H'
    ],
    'Omicron22B': [
        'G339D', 'S371F', 'S373P', 'S375F', 'T376A', 
        'D405N', 'K417N', 'N440K', 'L452R', 'S477N', 
        'T478K', 'E484A', 'F486V', 'Q498R', 'N501Y', 'Y505H'
    ],
    'Omicron22E': [
        'G339D', 'R346T', 'S371F', 'S373P', 'S375F', 
        'T376A', 'D405N', 'R408S', 'K417N', 'N440K', 
        'K444T', 'L452R', 'N460K', 'S477N', 'T478K', 
        'E484A', 'F486V', 'Q498R', 'N501Y', 'Y505H'
    ]
}

def mutation_to_residue(mutation):
    return int(''.join(filter(str.isdigit, mutation)))

def plot_mutations(ax, variant_mutations, color='gray', 
                  linestyle='--', alpha=0.7):
    for mutation in variant_mutations:
        residue = mutation_to_residue(mutation)
        ax.axvline(x=residue, color=color, linestyle=linestyle, 
                  alpha=alpha)

# Leer y procesar datos
results = {name: modify_xvg(read_xvg(Path(file)), shift) 
          for name, (file, shift) in files.items()}

# Calcular promedio para cada variante y encontrar residuos 
# por encima de la media
variant_means = {name: df['RMSF'].mean() 
                for name, df in results.items()}
above_average_residues = {}
for name, df in results.items():
    above_average_mask = df['RMSF'] > variant_means[name]
    above_average_residues[name] = set(
        df['Residue'][above_average_mask]
    )

# Encontrar residuos que están por encima de la media en todas 
# las variantes
consistently_high_residues = set.intersection(
    *above_average_residues.values()
)
dynamic_conserved_residues = sorted(consistently_high_residues)

# Graficar
fig, axes = plt.subplots(
    nrows=len(results), 
    ncols=1, 
    figsize=(10, 1.5*len(results)), 
    sharex=True
)
x_lim = (317, 541)
y_lim = (0, 0.68)
colors = ['black', 'blue', 'green', 'yellow', 'orange', 
          'lightcoral', 'coral']

for i, (name, df) in enumerate(results.items()):
    ax = axes[i]
    color = colors[i % len(colors)]
    ax.plot(df['Residue'], df['RMSF'], color=color)
    ax.fill_between(df['Residue'], df['RMSF'], color=color, 
                   alpha=0.5)
    
    for (start, end), highlight_color in highlight_regions.items():
        ax.axvspan(start, end, color=highlight_color, alpha=0.3)
    
    high_residue_mask = df['Residue'].isin(
        dynamic_conserved_residues
    )
    ax.scatter(
        df['Residue'][high_residue_mask],
        df['RMSF'][high_residue_mask],
        color='red', 
        s=20, 
        zorder=5
    )
    
    # Graficar la posición de las mutaciones en las variantes
    plot_mutations(ax, mutations[name])
    
    ax.set_ylim(y_lim)
    ax.set_ylabel('RMSF (nm)', fontsize=30)
    if i != len(results) // 2:
        ax.set_ylabel('')
    ax.tick_params(axis='y', labelsize=13)
    ax.text(
        1.02, 0.5, 
        name, 
        transform=ax.transAxes,
        fontsize=18, 
        fontweight='bold',
        verticalalignment='center',
        horizontalalignment='left'
    )

axes[-1].set_xlabel('Residuos', fontsize=30)
axes[-1].set_xlim(x_lim)
axes[-1].tick_params(axis='x', labelsize=12)

for ax in axes[:-1]:
    ax.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False
    )

plt.suptitle('', fontsize=18)
plt.tight_layout()
plt.subplots_adjust(right=0.85)
plt.show()

print("\nResiduos por encima del promedio en todas las variantes:", 
      dynamic_conserved_residues)
