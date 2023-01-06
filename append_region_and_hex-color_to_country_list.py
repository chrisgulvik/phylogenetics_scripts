#!/usr/bin/env python3

import sys

import pandas as pd
    # use openpyxl-3.0.9 for XLSX



report_first_region = True



# Read in the dataframe
xlsx = pd.ExcelFile('CIA-world-factbook-countries-by-region.xlsx')

# Store and print all sheet names for clarity
tabs = xlsx.sheet_names
for idx, tab in enumerate(tabs):
    sys.stderr.write('Sheet {}: \'{}\'\n'.format(idx, tab))

# Store and print all column names within the specific sheet
x1 = xlsx.parse('retrieved on 30-DEC-2022')
columns = x1.columns
for idx, col in enumerate(columns):
    sys.stderr.write('Column {}: \'{}\'\n'.format(idx, col))

# Store just the countries and regions as a dictionary
x1 = x1[['Country or Territory', 'Region']]
country_region_dict = x1.set_index('Country or Territory')['Region'].to_dict()

# Read in hex colors that match up with specific regions
df = pd.read_csv('hex-colors.tree-data.by-region.tsv', '\t')
colors_for_regions_dict = df.set_index('Region').to_dict('index')

# Read in line listing to be searched for country names
with open(sys.argv[1], 'r') as ifh:
    for line in ifh:
        line = str(line).rstrip('\n')
        ln = line.replace('_', ' ')
        region = list(map(country_region_dict.get, filter(lambda x:x in ln, country_region_dict)))
        # Jay often includes geo_loc_name as "USA_ex_foreign" to denote the
        #  foreign country travel history but that gives matches to different
        #  regions, so only save the last country after _ex_
        if ' ex ' in ln:
            _, l = ln.rsplit(' ex ', 1)
            region = list(map(country_region_dict.get, filter(lambda x:x in l, country_region_dict)))
        if len(region) == 1:
            d = colors_for_regions_dict[region[0]]
            s = '{}\t{}'.format(d['Tree Colors Data'], d['HEX Color'])
            print('{}\t{}\t{}'.format(line, s, region[0]))
        elif len(region) == 0:
            d = colors_for_regions_dict['Unknown']
            s = '{}\t{}'.format(d['Tree Colors Data'], d['HEX Color'])
            print('{}\t{}\t{}'.format(line, s, 'Unknown'))
        elif len(region) > 1 and report_first_region:
            d = colors_for_regions_dict[region[0]]
            s = '{}\t{}'.format(d['Tree Colors Data'], d['HEX Color'])
            print('{}\t{}\t{}'.format(line, s, region[0]))
        else:
            sys.stderr.write('ERROR: expected 0 or 1 country match to {} but found matches {}'.format(line, region))
            sys.exit(1)
