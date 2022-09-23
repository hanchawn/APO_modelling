# colours to use for each source/sim - taken from colorblind friendly palette below
colors = {'ecco': '#377eb8', 'jena': '#4daf4a', 'nemo': '#984ea3', 'ff': '#ff7f00', 'obs': '#f781bf'}
colors_sources = {f'ocean_{sim.lower()}' if sim.lower() in ['ecco', 'jena', 'nemo'] else sim.lower(): col
                  for sim, col in colors.items()}

# linestyles to use for each component
ls = {'o2': '-', 'co2': '--', 'n2': '-.', 'obs': 'none'}

# names of the models
model_names = {'ecco': 'ECCO', 'jena': 'Jena', 'nemo': 'NEMO-ERSEM'}

# adjustment to APO baseline
adjust = {2014: {1:0, 2:0, 3:15, 4:0, 5:0, 6:10, 7:5, 8:10, 9:0, 10:0, 11:10, 12:10},
          2015: {1:10, 2:0, 3:15, 4:10, 5:10, 6:10, 7:5, 8:10, 9:0, 10:0, 11:10, 12:5},
          2018: {1:0, 2:15, 3:40, 4:50, 5:50, 6:20, 7:0, 8:0, 9:15, 10:20, 11:20, 12:15},
          2020: {1:15, 2:30, 3:25, 4:40, 5:25, 6:0, 7:-5, 8:-10, 9:-5, 10:0, 11:0, 12:30}}

# colorblind friendly palette: https://gist.github.com/thriveth/8560036
# blue, orange, green,
# pink, brown, purple, 
# grey, red, yellow
CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

month_names = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
               7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}
seasons = {'Spring': [3, 4, 5], 'Summer': [6, 7, 8], 'Autumn': [9, 10, 11], 'Winter': [12, 1, 2]}
sitecodes = {'WAO': 'Weybourne', 'HFD': 'Heathfield', 'RGL': 'Ridge Hill'}
