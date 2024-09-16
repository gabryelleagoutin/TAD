#!/usr/bin/env python

import json
import pandas as pd
import argparse

__author__ = 'Gabryelle Agoutin - INRAE'
__copyright__ = 'Copyright (C) 2024 INRAE'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'gabryelle.agoutin@inrae.fr'
__status__ = 'prod'


def read_primers_from_table(filename):
    df = pd.read_csv(filename, delimiter='\t')
    return df

def generate_highcharts_script_from_primers(primers, chart_id_prefix, legend_prefix, letter_mapping, title):
    primer_colors = ['#277da1', '#577590', '#4d908e', '#43aa8b', '#90be6d', '#f9c74f', '#f9844a', '#f8961e', '#f3722c', '#f94144']

    series_data = []
    for i, primer in enumerate(primers, start=1):
        data = []
        for j, char in enumerate(primer, start=1):
            if char in letter_mapping:
                data.append({
                    'category': str(j),  #  Keep the position as a string
                    'y': letter_mapping[char]['y'],
                    'letter': letter_mapping[char]['letter']
                })
                
        series_data.append({
            'name': f'{legend_prefix} {i}',
            'id': f'{chart_id_prefix}_{i}',
            'color': primer_colors[i % len(primer_colors)],
            'data': data,
            'dataLabels': {
                'enabled': True,
                'color': '#000000',
                'format': '{point.letter}',
                'style': {
                    'fontSize': '10px',
                    'fontFamily': 'Verdana, sans-serif'
                }
            }
        })

    highcharts_script = f"""
Highcharts.chart('container{chart_id_prefix}', {{
    chart: {{
        type: 'line'
    }},
    title: {{
        text: '{title}'
    }},
    xAxis: {{
        categories: {json.dumps(list(range(1, len(primers[0]) + 2)))}, 
        labels: {{
            autoRotation: [-45, -90],
            step: 1,
            style: {{
                fontSize: '13px',
                fontFamily: 'Verdana, sans-serif'
            }}
        }},
        
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

def generate_xrange_chart_script(og_id, alignment_size, primers, primer_colors):
    series_data = []
    primer_color_map = {}

    for i, primer in enumerate(primers):
        color = primer_colors[i % len(primer_colors)]
        primer_color_map[primer['Index']] = color

        pos_A = primer['Position_A']
        primer_size_A = primer['Primer_Size_A']
        end_pos_A = pos_A + primer_size_A

        pos_B = primer['Position_B']
        primer_size_B = primer['Primer_Size_B']
        end_pos_B = pos_B + primer_size_B

        # Check whether a series for this primer index already exists
        serie_exists = False
        for serie in series_data:
            if serie['name'] == f'Primer {primer["Index"]}':
                serie['data'].append({
                    'x': pos_A,
                    'x2': end_pos_A,
                    'y': 0,
                    'color': color
                })
                serie['data'].append({
                    'x': pos_B,
                    'x2': end_pos_B,
                    'y': 1,
                    'color': color
                })
                serie_exists = True
                break

        if not serie_exists:
            series_data.append({
                'name': f'Primer {primer["Index"]}',
                'data': [{
                    'x': pos_A,
                    'x2': end_pos_A,
                    'y': 0,
                    'color': color
                }, {
                    'x': pos_B,
                    'x2': end_pos_B,
                    'y': 1,
                    'color': color
                }],
                'color': color,  
                'dataLabels': {
                    'enabled': True,
                    'color': '#000000',
                    'format': '{point.name}'
                }
            })

    highcharts_script = f"""
Highcharts.chart('container_{og_id}', {{
    chart: {{
        type: 'xrange'
    }},
    title: {{
        text: 'Primers position one the alignment for {og_id}'
    }},
    xAxis: {{
        title: {{
            text: 'Alignment Size'
        }},
        min: 0,
        max: {alignment_size}
    }},
    yAxis: {{
        categories: ['Forward Primer', 'Reverse Primer'],
        title: {{
            text: null
        }},
        reversed: true
    }},
    legend: {{
        enabled: true,
        itemStyle: {{
            color: '#000000'
        }},
        itemHoverStyle: {{
            color: '#FF0000'
        }},
        itemHiddenStyle: {{
            color: '#CCCCCC'
        }},
        symbolWidth: 30,
        symbolHeight: 10
    }},
    series: {json.dumps(series_data, indent=4)}
}});
"""
    return highcharts_script




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Generates an HTML file with Highcharts visualisations from a tabular primers file.""",
    epilog="Example:  python primer_metrics_visualization_prod.py -i sorted_results.tsv -o metrics_visualization.html")
    parser.add_argument('-i', '--input', required=True, help="Path to input file. This is the array of primers with selected pairs(sorted_results.tsv).")
    parser.add_argument('-o', '--output', required=True, help="output HTML name.")
    
    args = parser.parse_args()

    # df
    df = read_primers_from_table(args.input)

    # Highcharts Primer_A Reverse_Complement_B
    primers_A = df['Primer_A'].tolist()
    primers_B_reverse = df['Reverse_Complement_B'].tolist()

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

    highcharts_script_A = generate_highcharts_script_from_primers(primers_A, 'A', 'Primer', letter_mapping, 'Nucleotide Base Degeneracy Level')
    highcharts_script_B_reverse = generate_highcharts_script_from_primers(primers_B_reverse, 'B', 'Reverse_Primer', letter_mapping, 'Nucleotide Base Degeneracy Level')

    # min GC
    letter_mapping_gc = {
        'A': {'y': 1, 'letter': 'A'},
        'C': {'y': 2, 'letter': 'C'},
        'G': {'y': 2, 'letter': 'G'},
        'T': {'y': 1, 'letter': 'T'},
        'R': {'y': 1, 'letter': 'R'},
        'Y': {'y': 1, 'letter': 'Y'},
        'S': {'y': 2, 'letter': 'S'},
        'W': {'y': 1, 'letter': 'W'},
        'K': {'y': 1, 'letter': 'K'},
        'M': {'y': 1, 'letter': 'M'},
        'B': {'y': 1, 'letter': 'B'},
        'D': {'y': 1, 'letter': 'D'},
        'H': {'y': 1, 'letter': 'H'},
        'V': {'y': 1, 'letter': 'V'},
        'N': {'y': 1, 'letter': 'N'}
    }

    highcharts_script_A_gc = generate_highcharts_script_from_primers(primers_A, 'A_GC', 'Primer_GC', letter_mapping_gc, 'Nucleotide Minimal GC Content')
    highcharts_script_B_reverse_gc = generate_highcharts_script_from_primers(primers_B_reverse, 'B_GC', 'Reverse_Primer_GC', letter_mapping_gc, 'Nucleotide Minimal GC Content')

    # max GC
    letter_mapping_gc_max = {
        'A': {'y': 1, 'letter': 'A'},
        'C': {'y': 2, 'letter': 'C'},
        'G': {'y': 2, 'letter': 'G'},
        'T': {'y': 1, 'letter': 'T'},
        'R': {'y': 2, 'letter': 'R'},
        'Y': {'y': 2, 'letter': 'Y'},
        'S': {'y': 2, 'letter': 'S'},
        'W': {'y': 1, 'letter': 'W'},
        'K': {'y': 2, 'letter': 'K'},
        'M': {'y': 1, 'letter': 'M'},
        'B': {'y': 2, 'letter': 'B'},
        'D': {'y': 2, 'letter': 'D'},
        'H': {'y': 2, 'letter': 'H'},
        'V': {'y': 2, 'letter': 'V'},
        'N': {'y': 1, 'letter': 'N'}
    }

    highcharts_script_A_gc_max = generate_highcharts_script_from_primers(primers_A, 'A_GC_Max', 'Primer_GC_Max', letter_mapping_gc_max, 'Nucleotide Maximal GC Content')
    highcharts_script_B_reverse_gc_max = generate_highcharts_script_from_primers(primers_B_reverse, 'B_GC_Max', 'Reverse_Primer_GC_Max', letter_mapping_gc_max, 'Nucleotide Maximal GC Content')

    # x-range
    og_ids = df['OG_ID'].unique()
    primer_colors = ['#277da1', '#577590', '#4d908e', '#43aa8b', '#90be6d', '#f9c74f', '#f9844a', '#f8961e', '#f3722c', '#f94144']
    highcharts_scripts_xrange = []

    for og_id in og_ids:
        og_data = df[df['OG_ID'] == og_id]
        alignment_size = og_data['Alignement_size'].iloc[0]
        primers = []
        for i, row in og_data.iterrows():
            primers.append({
                'Index': i + 1,
                'Position_A': row['Position_A'],
                'Primer_Size_A': row['Primer_Size_A'],
                'Position_B': row['Position_B'],
                'Primer_Size_B': row['Primer_Size_B']
            })
        script = generate_xrange_chart_script(og_id, alignment_size, primers, primer_colors)
        highcharts_scripts_xrange.append(script)

    # Template HTML
    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Visualisation of primers metrics</title>
        <script src="https://code.highcharts.com/highcharts.js"></script>
        <script src="https://code.highcharts.com/modules/xrange.js"></script>
        <script src="https://code.highcharts.com/modules/exporting.js"></script>
        <script src="https://code.highcharts.com/modules/export-data.js"></script>
        <script src="https://code.highcharts.com/modules/accessibility.js"></script>
        <style>
            body {{
                text-align: center; 
            }}
            h1, h2 {{
                text-align: center; 
            }}
            .highcharts-figure {{
                min-width: 310px;
                max-width: 800px;
                margin: 1em auto; 
            }}
            #container {{
                height: 300px;
            }}
        </style>
    </head>
    <body>
    <h1>Visualisation of primers metrics</h1>
    <!-- Graphiques linéaires -->
    <h2>1. Degeneracy</h2>
    <figure class="highcharts-figure">
        <div id="containerA"></div>
        <p class="highcharts-description">Forward</p>
    </figure>
    <figure class="highcharts-figure">
        <div id="containerB"></div>
        <p class="highcharts-description">Reverse</p>
    </figure>

    <h2>2. GC Content</h2>
    <figure class="highcharts-figure">
        <div id="containerA_GC"></div>
        <p class="highcharts-description">Forward</p>
    </figure>
    <figure class="highcharts-figure">
        <div id="containerB_GC"></div>
        <p class="highcharts-description">Reverse</p>
    </figure>

    <figure class="highcharts-figure">
        <div id="containerA_GC_Max"></div>
        <p class="highcharts-description">Forward</p>
    </figure>
    <figure class="highcharts-figure">
        <div id="containerB_GC_Max"></div>
        <p class="highcharts-description">Reverse</p>
    </figure>

    <!-- Graphiques X-range -->
    <h2>3. Primer Position</h2>
    {''.join(f'<figure class="highcharts-figure"><div id="container_{og_id}"></div><p class="highcharts-description">Graphique pour {og_id}</p></figure>' for og_id in og_ids)}

    <script>
        (function(H) {{
            H.Legend.prototype.getAllItems = function () {{
                var allItems = [];
                H.each(this.chart.series, function (series) {{
                    var seriesOptions = series && series.options;
                    if (series.type === 'xrange') {{
                        series.color = series.userOptions.color;
                    }}
                    if (series && H.pick(
                        seriesOptions.showInLegend,
                        !H.defined(seriesOptions.linkedTo) ? undefined : false, true
                    )) {{
                        allItems = allItems.concat(
                            series.legendItems ||
                            (
                                seriesOptions.legendType === 'point' ?
                                    series.data :
                                    series
                            )
                        );
                     }}
                }});
                H.fireEvent(this, 'afterGetAllItems', {{ allItems: allItems }});
                return allItems;
            }}
        }})(Highcharts);

        {highcharts_script_A}
        {highcharts_script_B_reverse}
        {highcharts_script_A_gc}
        {highcharts_script_B_reverse_gc}
        {highcharts_script_A_gc_max}
        {highcharts_script_B_reverse_gc_max}

        {''.join(script for script in highcharts_scripts_xrange)}
    </script>
    </body>
    </html>
    """


    with open(args.output, 'w') as file:
        file.write(html_template)

    print(f"{args.output} file generated successfully.")


