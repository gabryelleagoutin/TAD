import json
import pandas as pd

def read_primers_from_table(filename):
    # Charger le fichier tabulaire en utilisant Pandas
    df = pd.read_csv(filename, delimiter='\t')
    
    # Extraire les primers de la colonne 'Primer_A' et 'Reverse_Complement_B'
    primers_A = df['Primer_A'].tolist()
    primers_B_reverse = df['Reverse_Complement_B'].tolist()
    
    return primers_A, primers_B_reverse

def generate_highcharts_script_from_primers(primers, chart_id_prefix, legend_prefix):
    # Correspondance des lettres avec leur catégorie et représentation
    letter_mapping = {
        'A': {'y': 1, 'letter': 'A'},
        'C': {'y': 1, 'letter': 'C'},
        'G': {'y': 1, 'letter': 'G'},
        'T': {'y': 1, 'letter': 'T'},
        'R': {'y': 2, 'letter': 'R'},
        'Y': {'y': 2, 'letter': 'Y'},
        'S': {'y': 2, 'letter': 'S'},
        'W': {'y': 2, 'letter': 'W'},
        'K': {'y': 2, 'letter': 'K'},
        'M': {'y': 2, 'letter': 'M'},
        'B': {'y': 3, 'letter': 'B'},
        'D': {'y': 3, 'letter': 'D'},
        'H': {'y': 3, 'letter': 'H'},
        'V': {'y': 3, 'letter': 'V'},
        'N': {'y': 4, 'letter': 'N'}
    }

    # Liste des couleurs pour chaque primer
    primer_colors = ['#277da1', '#577590', '#4d908e', '#43aa8b', '#90be6d', '#f9c74f', '#f9844a', '#f8961e', '#f3722c', '#f94144']

    series_data = []
    for i, primer in enumerate(primers, start=1):
        data = []
        for j, char in enumerate(primer, start=1):
            if char in letter_mapping:
                data.append({
                    'category': str(j),
                    'y': letter_mapping[char]['y'],
                    'letter': letter_mapping[char]['letter']
                })
        series_data.append({
            'name': f'{legend_prefix} {i}',
            'id': f'{chart_id_prefix}_{i}',  # Utilisation d'un identifiant unique pour chaque série
            'color': primer_colors[i % len(primer_colors)],  # Assignation de couleur par primer
            'groupPadding': 0.2,
            'pointWidth': 5,
            'data': data,
            'dataLabels': {
                'enabled': False,
                'rotation': -90,
                'color': '#FFFFFF',
                'inside': True,
                'verticalAlign': 'top',
                'format': '{point.y:.1f}',
                'y': 10,
                'style': {
                    'fontSize': '10px',
                    'fontFamily': 'Verdana, sans-serif'
                }
            }
        })

    highcharts_script = f"""
Highcharts.chart('container{chart_id_prefix}', {{
    chart: {{
        type: 'column'
    }},
    title: {{
        text: 'Nucleotide Base Degeneracy Level'
    }},
    xAxis: {{
        type: 'category',
        labels: {{
            autoRotation: [-45, -90],
            style: {{
                fontSize: '13px',
                fontFamily: 'Verdana, sans-serif'
            }}
        }}
    }},
    yAxis: {{
        min: 0,
        title: {{
            text: 'Degeneracy Level'
        }}
    }},
    legend: {{
        enabled: true
    }},
    tooltip: {{
        formatter: function () {{
            return '<b>' + this.series.name + '</b><br/>' +
                   'Position ' + this.point.category + ': ' +
                   this.point.options.letter + '<br/>' +
                   'Valeur: ' + this.point.y;
        }}
    }},
    series: {json.dumps(series_data, indent=4)}
}});
"""
    return highcharts_script

# Lecture des primers depuis les colonnes 'Primer_A' et 'Reverse_Complement_B'
primers_A, primers_B_reverse = read_primers_from_table('sorted_results.tsv')

# Génération du script Highcharts pour Primer_A
highcharts_script_A = generate_highcharts_script_from_primers(primers_A, 'A', 'Primer')

# Génération du script Highcharts pour Reverse_Complement_B
highcharts_script_B_reverse = generate_highcharts_script_from_primers(primers_B_reverse, 'B', 'Reverse_Primer')

# Template HTML
html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Highcharts Example</title>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
    <script src="https://code.highcharts.com/modules/export-data.js"></script>
    <script src="https://code.highcharts.com/modules/accessibility.js"></script>
    <style>
        .highcharts-figure {{
            min-width: 310px;
            max-width: 800px;
            margin: 1em auto;
        }}
    </style>
</head>
<body>

<figure class="highcharts-figure">
    <div id="containerA"></div>
    <p class="highcharts-description">Graphique pour Primer_A</p>
</figure>

<figure class="highcharts-figure">
    <div id="containerB"></div>
    <p class="highcharts-description">Graphique pour Reverse_Complement_B</p>
</figure>

<script>
    {highcharts_script_A}
    {highcharts_script_B_reverse}
    
    Highcharts.chart('container', {{
        plotOptions: {{
            column: {{
                borderWidth: 0  // Supprime les bordures autour des colonnes
            }}
        }}
    }});    
</script>

</body>
</html>
"""

# Écriture dans un fichier HTML
with open('highcharts_example.html', 'w') as file:
    file.write(html_template)

print("HTML file generated successfully.")
