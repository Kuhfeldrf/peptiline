import matplotlib.pyplot as plt

# Single color options with category headers
single_color_options = [
    '--- SINGLE COLORS OPTIONS ---',
    'red', 'green', 'blue', 'yellow', 'purple', 'orange', 'cyan', 
    'magenta', 'pink', 'brown', 'black', 'white', 'gray', 'darkblue',
    'darkgreen', 'darkred', 'darkorange', 'darkviolet', 'lightblue',
    'lightgreen', 'lightcoral', 'gold', 'silver', 'teal', 'navy', 'maroon',
    'olive', 'lime', 'aqua', 'indigo', 'violet', 'turquoise', 'coral',
    'crimson', 'salmon', 'sienna', 'tan', 'khaki', 'plum', 'orchid'
]

# Discrete colormaps with category header
valid_discrete_cmaps = [
    '--- QUALITATIVE PALETTES ---',
    'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 
    'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c'
]

# Gradient colormaps organized with category headers
valid_sequential_cmaps = [
    '--- SEQUENTIAL PALETTES ---',
    'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'OrRd', 
    'Oranges', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'Reds', 'YlGn', 
    'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 
    'cool', 'copper', 'cubehelix', 'gist_earth', 'gist_gray', 'gist_heat', 
    'gist_yarg', 'gray', 'hot', 'pink', 'summer', 'winter', 'grey', 
    'gist_grey', 'gist_yerg', 'Grays'
]

valid_diverging_cmaps = [
    '--- DIVERGING PALETTES ---',
    'BrBG', 'PRGn', 'PiYG', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 
    'Spectral', 'bwr', 'coolwarm', 'seismic'
]

valid_cyclical_cmaps = [
    '--- CYCLICAL PALETTES ---',
    'jet', 'ocean', 'rainbow', 'terrain', 'brg', 'CMRmap', 
    'flag', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gnuplot', 
    'gnuplot2', 'nipy_spectral', 'prism', 'Wistia'
]

# Dictionary for reference (not for direct dropdown use)
valid_gradient_cmaps = {
    'sequential': [c for c in valid_sequential_cmaps if not c.startswith('---')],
    'diverging': [c for c in valid_diverging_cmaps if not c.startswith('---')],
    'cyclical': [c for c in valid_cyclical_cmaps if not c.startswith('---')]
}

# Combined list for dropdowns with all categories
all_gradient_cmaps = [
    '--- DEFAULT PALETTE(RdYlGn_r) ---',
    'RdYlGn_r'
]
all_gradient_cmaps.extend(valid_sequential_cmaps)
all_gradient_cmaps.extend(valid_diverging_cmaps)
all_gradient_cmaps.extend(valid_cyclical_cmaps)

# Add reversed versions with their own category
all_gradient_cmaps.append('--- REVERSED PALETTES ---')
for category in valid_gradient_cmaps.values():
    for cmap in category:
        all_gradient_cmaps.append(cmap + '_r')



# Species translation list for proteomic and peptidomic analysis
SPEC_TRANSLATE_LIST = [
    # Primary model organisms
    ['Human', 'homo sapiens', 'human', 'h. sapiens', '9606'],
    ['Mouse', 'mus musculus', 'murine', 'm. musculus', '10090'],
    ['Rat', 'rattus norvegicus', 'r. norvegicus', '10116'],
    
    # Livestock/agriculture 
    ['Bovine', 'bos taurus', 'cow', 'cattle', 'b. taurus', '9913'],
    ['Pig', 'sus scrofa', 'porcine', 's. scrofa', 'swine', '9823'],
    ['Sheep', 'ovis aries', 'o. aries', '9940'],
    ['Goat', 'capra hircus', 'c. hircus', '9925'],
    ['Chicken', 'gallus gallus', 'g. gallus', '9031'],
    ['Turkey', 'meleagris gallopavo', 'm. gallopavo', '9103'],
    ['Duck', 'anas platyrhynchos', 'a. platyrhynchos', '8839'],
    ['Quail', 'coturnix japonica', 'c. japonica', '93934'],
    
    # Common lab organisms
    ['Fruit fly', 'drosophila melanogaster', 'd. melanogaster', '7227'],
    ['Zebrafish', 'danio rerio', 'd. rerio', '7955'],
    ['Roundworm', 'caenorhabditis elegans', 'c. elegans', '6239'],
    ['Yeast', 'saccharomyces cerevisiae', 's. cerevisiae', '4932'],
    ['E. coli', 'escherichia coli', 'e. coli', '562'],
    ['Arabidopsis', 'arabidopsis thaliana', 'a. thaliana', 'thale cress', '3702'],
    ['Fission yeast', 'schizosaccharomyces pombe', 's. pombe', '4896'],
    ['Xenopus', 'xenopus laevis', 'x. laevis', 'african clawed frog', '8355'],
    
    # Primates
    ['Human', 'homo sapiens', 'human', 'h. sapiens', '9606'],
    ['Chimpanzee', 'pan troglodytes', 'p. troglodytes', '9598'],
    ['Bonobo', 'pan paniscus', 'p. paniscus', '9597'],
    ['Gorilla', 'gorilla gorilla', 'g. gorilla', '9593'],
    ['Orangutan', 'pongo abelii', 'p. abelii', '9601'],
    ['Rhesus macaque', 'macaca mulatta', 'm. mulatta', '9544'],
    ['Marmoset', 'callithrix jacchus', 'c. jacchus', '9483'],
    ['Baboon', 'papio anubis', 'p. anubis', '9555'],
    
    # Other mammals
    ['Horse', 'equus caballus', 'e. caballus', 'equine', '9796'],
    ['Rabbit', 'oryctolagus cuniculus', 'o. cuniculus', '9986'],
    ['Dog', 'canis lupus familiaris', 'canis familiaris', 'c. familiaris', '9615'],
    ['Cat', 'felis catus', 'f. catus', '9685'],
    ['Guinea pig', 'cavia porcellus', 'c. porcellus', '10141'],
    ['Hamster', 'mesocricetus auratus', 'm. auratus', '10036'],
    ['Ferret', 'mustela putorius furo', 'm. putorius', '9669'],
    ['Hedgehog', 'erinaceus europaeus', 'e. europaeus', '9365'],
    
    # Marine organisms
    ['Zebrafish', 'danio rerio', 'd. rerio', '7955'],
    ['Atlantic salmon', 'salmo salar', 's. salar', '8030'],
    ['Rainbow trout', 'oncorhynchus mykiss', 'o. mykiss', '8022'],
    ['Tilapia', 'oreochromis niloticus', 'o. niloticus', '8128'],
    ['Cod', 'gadus morhua', 'g. morhua', '8049'],
    ['Dolphin', 'tursiops truncatus', 't. truncatus', 'bottlenose dolphin', '9739'],
    ['Killer whale', 'orcinus orca', 'o. orca', '9733'],
    ['Sea urchin', 'strongylocentrotus purpuratus', 's. purpuratus', '7668'],
    
    # Microorganisms
    ['E. coli', 'escherichia coli', 'e. coli', '562'],
    ['B. subtilis', 'bacillus subtilis', 'b. subtilis', '1423'],
    ['S. aureus', 'staphylococcus aureus', 's. aureus', '1280'],
    ['P. aeruginosa', 'pseudomonas aeruginosa', 'p. aeruginosa', '287'],
    ['M. tuberculosis', 'mycobacterium tuberculosis', 'm. tuberculosis', '1773'],
    ['L. monocytogenes', 'listeria monocytogenes', 'l. monocytogenes', '1639'],
    ['S. pneumoniae', 'streptococcus pneumoniae', 's. pneumoniae', '1313'],
    ['S. cerevisiae', 'saccharomyces cerevisiae', 's. cerevisiae', 'baker\'s yeast', '4932'],
    ['S. pombe', 'schizosaccharomyces pombe', 's. pombe', 'fission yeast', '4896'],
    ['C. albicans', 'candida albicans', 'c. albicans', '5476'],
    
    # Plants
    ['Arabidopsis', 'arabidopsis thaliana', 'a. thaliana', 'thale cress', '3702'],
    ['Rice', 'oryza sativa', 'o. sativa', '4530'],
    ['Maize', 'zea mays', 'z. mays', 'corn', '4577'],
    ['Wheat', 'triticum aestivum', 't. aestivum', '4565'],
    ['Soybean', 'glycine max', 'g. max', '3847'],
    ['Tomato', 'solanum lycopersicum', 's. lycopersicum', '4081'],
    ['Potato', 'solanum tuberosum', 's. tuberosum', '4113'],
    ['Tobacco', 'nicotiana tabacum', 'n. tabacum', '4097'],
    
    # Insects
    ['Fruit fly', 'drosophila melanogaster', 'd. melanogaster', '7227'],
    ['Honeybee', 'apis mellifera', 'a. mellifera', '7460'],
    ['Mosquito', 'anopheles gambiae', 'a. gambiae', '7165'],
    ['Silkworm', 'bombyx mori', 'b. mori', '7091'],
    ['Beetle', 'tribolium castaneum', 't. castaneum', 'red flour beetle', '7070'],
    
    # Additional species from original list
    ['Buffalo', 'bubalus bubalis', 'water buffalo', 'b. bubalis', '89462'],
    ['Donkey', 'equus asinus', 'e. asinus', '9793'],
    ['Camel', 'camelus dromedarius', 'dromedary camel', 'c. dromedarius', '9838'],
    ['Yak', 'bos mutus', 'b. mutus', '72004'],
    ['Deer', 'cervus elaphus', 'red deer', 'c. elaphus', '9860'],
    ['Reindeer', 'rangifer tarandus', 'r. tarandus', '9870'],
    ['American bison', 'bison bison', 'b. bison', '9901'],
    ['Elephant', 'elephas maximus', 'asian elephant', 'e. maximus', '9783']
]