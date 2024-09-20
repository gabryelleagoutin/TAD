import pandas as pd
import json
import argparse

def read_primers_from_table(filename):
    df = pd.read_csv(filename, delimiter='\t')
    return df

def generate_split_yaxis_highcharts_script(primers_A, primers_B, chart_id_prefix, title):
    letter_mapping = {
        'A': 1, 'C': 1, 'G': 1, 'T': 1,
        'R': 2, 'Y': 2, 'S': 2, 'W': 2, 'K': 2, 'M': 2,
        'B': 3, 'D': 3, 'H': 3, 'V': 3, 'N': 4
    }

    series_A = []
    series_B = []
    for j, char in enumerate(primers_A, start=1):
        series_A.append({
            'category': str(j),
            'y': letter_mapping.get(char, 1),
            'letter': char
        })

    for j, char in enumerate(primers_B, start=1):
        series_B.append({
            'category': str(j),
            'y': letter_mapping.get(char, 1),
            'letter': char
        })

    split_highcharts_script = f"""
Highcharts.chart('container{chart_id_prefix}', {{
    chart: {{
        type: 'line',
        spacingBottom: 50
    }},
    title: {{
        text: '{title}'
    }},
    xAxis: {{
        categories: {json.dumps(list(range(1, len(primers_A) + 1)))}
    }},
    yAxis: [{{
        title: {{
            text: 'Degeneracy Level',
            align: 'middle',
            rotation: 270,
            margin: 40
        }},
        min: 0,
        max: 4,
        tickPositions: [0, 1, 2, 3, 4],
        labels: {{
            align: 'left',
            x: -10,
            y: 10,
            style: {{
                fontSize: '10px'
            }}
        }},
        height: '45%',
        offset: 0
    }}, {{
        title: {{
            text: 'Degeneracy Level',
            align: 'middle',
            rotation: 270,
            margin: 40
        }},
        min: 0,
        max: 4,
        tickPositions: [0, 1, 2, 3, 4],
        labels: {{
            align: 'left',
            x: -10,
            y: 10,
            style: {{
                fontSize: '10px'
            }}
        }},
        top: '55%',
        height: '45%',
        offset: 0
    }}],
    plotOptions: {{
        series: {{
            dataLabels: {{
                enabled: true,
                format: '{{point.letter}}',
                style: {{
                    fontSize: '10px'
                }}
            }}
        }}
    }},
    series: [{{
        name: 'Forward Primer',
        data: {json.dumps(series_A)},
        yAxis: 0
    }}, {{
        name: 'Reverse Primer',
        data: {json.dumps(series_B)},
        yAxis: 1
    }}]
}});
"""
    return split_highcharts_script

def generate_split_gc_content_script(primers_A, primers_B, chart_id_prefix, title):
    gc_mapping = {'A': 0, 'C': 1, 'G': 1, 'T': 0, 'R': 0.5, 'Y': 0.5, 'S': 1, 'W': 0, 'K': 0.5, 'M': 0.5, 'B': 0.667, 'D': 0.333, 'H': 0.333, 'V': 0.667, 'N': 0.5}

    series_A = []
    series_B = []

    for j, char in enumerate(primers_A, start=1):
        series_A.append({
            'category': str(j),
            'y': gc_mapping.get(char, 0),
            'letter': char
        })

    for j, char in enumerate(primers_B, start=1):
        series_B.append({
            'category': str(j),
            'y': gc_mapping.get(char, 0),
            'letter': char
        })

    split_gc_content_script = f"""
Highcharts.chart('container{chart_id_prefix}', {{
    chart: {{
        type: 'line',
        spacingTop: 50
    }},
    title: {{
        text: '{title}'
    }},
    xAxis: {{
        categories: {json.dumps(list(range(1, len(primers_A) + 1)))}
    }},
    yAxis: [{{
        title: {{
            text: 'GC rate',
            align: 'middle',
            rotation: 270,
            margin: 40
        }},
        min: 0.1,
        max: 1.5,
        tickPositions: [0, 0.5, 1, 1.5],
        labels: {{
            align: 'left',
            x: -10,
            y: 10,
            style: {{
                fontSize: '10px'
            }}
        }},
        height: '45%',
        offset: 0
    }}, {{
        title: {{
            text: 'GC rate',
            align: 'middle',
            rotation: 270,
            margin: 40
        }},
        min: 0.1,
        max: 1.5,
        tickPositions: [0, 0.5, 1, 1.5],
        labels: {{
            align: 'left',
            x: -10,
            y: 10,
            style: {{
                fontSize: '10px'
            }}
        }},
        top: '55%',
        height: '45%',
        offset: 0
    }}],
    plotOptions: {{
        series: {{
            dataLabels: {{
                enabled: true,
                format: '{{point.letter}}',
                style: {{
                    fontSize: '10px'
                }}
            }}
        }}
    }},
    series: [{{
        name: 'Forward Primer',
        data: {json.dumps(series_A)},
        yAxis: 0
    }}, {{
        name: 'Reverse Primer',
        data: {json.dumps(series_B)},
        yAxis: 1
    }}]
}});
"""
    return split_gc_content_script

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

        series_data.append({
            'name': f'Primer pair {primer["Index"]}',
            'data': [
                {'x': pos_A, 'x2': end_pos_A, 'y': 0, 'color': color},
                {'x': pos_B, 'x2': end_pos_B, 'y': 1, 'color': color}
            ],
            'color': color
        })

    highcharts_script = f"""
Highcharts.chart('container_xrange_{og_id}', {{
    chart: {{
        type: 'xrange'
    }},
    title: {{
        text: 'Primer Position for {og_id}'
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
    series: {json.dumps(series_data, indent=4)}
}});
"""
    return highcharts_script

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Generates an HTML file with Highcharts visualizations from a tabular primers file.""",
                                     epilog="Example: python primer_metrics_visualization.py -i sorted_results.tsv -o visualization.html")

    parser.add_argument('-i', '--input', required=True, help="Path to input file. This is the array of primers with selected pairs (sorted_results.tsv).")
    parser.add_argument('-o', '--output', required=True, help="Output HTML file name.")

    args = parser.parse_args()

    df = read_primers_from_table(args.input)

    # Initialize variables to collect all chart scripts and chart divs
    all_chart_scripts = ""
    all_content_blocks = ""

    # Generate the existing line and GC content charts
    for index, row in df.iterrows():
        primer_data = row.to_dict()

        chart_id_prefix = f"Combined_{index}"
        gc_chart_id_prefix = f"GC_Combined_{index}"

        primers_A = row['Primer_A']
        primers_B_reverse = row['Reverse_Complement_B']

        split_chart_script = generate_split_yaxis_highcharts_script(primers_A, primers_B_reverse, chart_id_prefix, 'Forward and Reverse Primer Degeneracy Level')
        split_gc_content_script = generate_split_gc_content_script(primers_A, primers_B_reverse, gc_chart_id_prefix, 'Forward and Reverse Primer GC Content')

        all_chart_scripts += split_chart_script + "\n" + split_gc_content_script + "\n"

        # Add a title for each pair of primers (Primer Pair 1, Primer Pair 2, etc.)
        pair_title = f"Primer pair {index + 1}"

        content_block = f"""
        <h2>{pair_title}</h2>  <!-- Adding dynamic title for each primer pair -->
        <div class="content-wrapper">
            <div class="table-chart-container">
                <table>
                    <tr>
                        <th colspan="2">OG ID</th>
                        <td colspan="2">{primer_data.get('OG_ID', '')}</td>
                    </tr>
                    <tr>
                        <th colspan="2">Gene Name</th>
                        <td colspan="2">{primer_data.get('GeneName', '')}</td>
                    </tr>
                    <tr>
                        <th colspan="2">Potential Amplicon Size</th>
                        <td colspan="2">{primer_data.get('potential_amplicon_size', '')}</td>
                    </tr>
                    <tr>
                        <th colspan="2">TxMk Score</th>
                        <td colspan="2">{primer_data.get('Total_score', '')}</td>
                    </tr>
                    <tr>
                        <th colspan="2">Min Size Amplicon</th>
                        <td colspan="2">{primer_data.get('min_size_amplicon', '')}</td>
                    </tr>
                    <tr>
                        <th colspan="2">Max Size Amplicon</th>
                        <td colspan="2">{primer_data.get('max_size_amplicon', '')}</td>
                    </tr>
                    <tr>
                        <th colspan="2" style="text-align: center;">Forward Primer</th>
                        <th colspan="2" style="text-align: center;">Reverse Primer</th>
                    </tr>
                    <tr>
                        <th>sequence</th><td>{primer_data.get('Primer_A', '')}</td>
                        <th>sequence</th><td>{primer_data.get('Primer_B', '')}</td>
                    </tr>
                    <tr>
                        <th>Position</th><td>{primer_data.get('Position_A', '')}</td>
                        <th>Position</th><td>{primer_data.get('Position_B', '')}</td>
                    </tr>
                    <tr>
                        <th>Number Matching</th><td>{primer_data.get('Number_matching_A', '')}</td>
                        <th>Number Matching</th><td>{primer_data.get('Number_matching_B', '')}</td>
                    </tr>
                    <tr>
                        <th>Percentage NM</th><td>{primer_data.get('Percentage_NM_A', '')}</td>
                        <th>Percentage NM</th><td>{primer_data.get('Percentage_NM_B', '')}</td>
                    </tr>
                    <tr>
                        <th>Degenerescence</th><td>{primer_data.get('Degenerescence_A', '')}</td>
                        <th>Degenerescence</th><td>{primer_data.get('Degenerescence_B', '')}</td>
                    </tr>
                    <tr>
                        <th>Tm Max (°C)</th><td>{primer_data.get('Tm_A_max', '')}</td>
                        <th>Tm Max (°C)</th><td>{primer_data.get('Tm_B_max', '')}</td>
                    </tr>
                    <tr>
                        <th>Tm Min (°C)</th><td>{primer_data.get('Tm_A_min', '')}</td>
                        <th>Tm Min (°C)</th><td>{primer_data.get('Tm_B_min', '')}</td>
                    </tr>
                    <tr>
                        <th>GC Percentage Fraction</th><td>{primer_data.get('GC_percentage_fraction_A', '')}</td>
                        <th>GC Percentage Fraction</th><td>{primer_data.get('GC_percentage_fraction_B', '')}</td>
                    </tr>
                    <tr>
                        <th>GC in Last 30%</th><td>{primer_data.get('GC_in_last_thirty_percent_A', '')}</td>
                        <th>GC in Last 30%</th><td>{primer_data.get('GC_in_last_thirty_percent_B', '')}</td>
                    </tr>
                    <tr>
                        <th>Ends with T </th><td>{primer_data.get('Ends_with_T_A', '')}</td>
                        <th>Ends with T </th><td>{primer_data.get('Ends_with_T_B', '')}</td>
                    </tr>
                    <tr>
                        <th>GC Clamp </th><td>{primer_data.get('GC_clamp_A', '')}</td>
                        <th>GC Clamp </th><td>{primer_data.get('GC_clamp2', '')}</td>
                    </tr>
                    <tr>
                        <th>Reverse Complement </th><td>-</td>
                        <th>Reverse Complement </th><td>{primer_data.get('Reverse_Complement_B', '')}</td>
                    </tr>
                    <tr>
                        <th>GC Clamp RC </th><td>-</td>
                        <th>GC Clamp RC </th><td>{primer_data.get('GC_clamp_RC_B', '')}</td>
                    </tr>
                </table>
                
                <div class="chart-container">
                    <div class="chart-column">
                        <div id="container{chart_id_prefix}"></div>
                    </div>
                    <div class="chart-column">
                        <div id="container{gc_chart_id_prefix}"></div>
                    </div>
                </div>
            </div>
        </div>
        """

        all_content_blocks += content_block

    # Add the X-range chart for Primer Position
    og_ids = df['OG_ID'].unique()
    primer_colors = ['#277da1', '#577590', '#4d908e', '#43aa8b', '#90be6d', '#f9c74f', '#f9844a', '#f8961e', '#f3722c', '#f94144']
    xrange_scripts = ""
    
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
        xrange_scripts += generate_xrange_chart_script(og_id, alignment_size, primers, primer_colors) + "\n"

    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Visualization of Primer Metrics</title>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/xrange.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
    <script src="https://code.highcharts.com/modules/export-data.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            text-align: center;
        }}
        .content-wrapper {{
            border: 2px solid black;
            padding: 20px;
            display: inline-block;
            margin-bottom: 50px;
            width: 90%;
        }}
        .table-chart-container {{
            display: flex;
        }}
        table {{
            width: 50%;  
            font-size: 12px;
            margin: 10px;
            border-collapse: collapse;
        }}
        table, th, td {{
            border: 1px solid black;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
        }}
        .chart-container {{
            display: flex;
            flex-direction: column;
            justify-content: space-between;
            width: 50%; 
            margin-left: 10px;
        }}
        .chart-column {{
            width: 100%;
            height: 400px;
        }}
    </style>
</head>
<body>
    <h1>Primer Visualization</h1>

    {all_content_blocks}

    <!-- X-range -->
    <h2>Primer Position</h2>
    {''.join(f'<div id="container_xrange_{og_id}"></div>' for og_id in og_ids)}

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

        {all_chart_scripts}
        {xrange_scripts}
    </script>
</body>
</html>
"""

    # Write to the output file
    with open(args.output, 'w') as f:
        f.write(html_template)

    print(f"HTML file generated: {args.output}")
