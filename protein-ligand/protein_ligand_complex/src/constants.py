# all punctuation
punctuation_regex  = r"""(\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"""

# tokenization regex (Schwaller)
molecule_regex = r"""(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"""

# filter out these common additives which occur in more than 75 complexes in the PDB
ubiquitous_ligands = [
    'PEG', 'ADP', 'FAD', 'NAD', 'ATP', 'MPD', 'NAP', 'GDP', 'MES',
    'GTP', 'FMN', 'HEC', 'TRS', 'CIT', 'PGE', 'ANP', 'SAH', 'NDP',
    'PG4', 'EPE', 'AMP', 'COA', 'MLI', 'FES', 'GNP', 'MRD', 'GSH',
    'FLC', 'AGS', 'NAI', 'SAM', 'PCW', '1PE', 'TLA', 'BOG', 'CYC',
    'UDP', 'PX4', 'NAG', 'IMP', 'POP', 'UMP', 'PLM', 'HEZ', 'TPP',
    'ACP', 'LDA', 'ACO', 'CLR', 'BGC', 'P6G', 'LMT', 'OGA', 'DTT',
    'POV', 'FBP', 'AKG', 'MLA', 'ADN', 'NHE', '7Q9', 'CMP', 'BTB',
    'PLP', 'CAC', 'SIN', 'C2E', '2AN', 'OCT', '17F', 'TAR', 'BTN',
    'XYP', 'MAN', '5GP', 'GAL', 'GLC', 'DTP', 'DGT', 'PEB', 'THP',
    'BEZ', 'CTP', 'GSP', 'HED', 'ADE', 'TYD', 'TTP', 'BNG', 'IHP',
    'FDA', 'PEP', 'ALF', 'APR', 'MTX', 'MLT', 'LU8', 'UTP', 'APC',
    'BLA', 'C8E', 'D10', 'CHT', 'BO2', '3BV', 'ORO', 'MPO', 'Y01',
    'OLC', 'B3P', 'G6P', 'PMP', 'D12', 'NDG', 'A3P', '78M', 'F6P',
    'U5P', 'PRP', 'UPG', 'THM', 'SFG', 'MYR', 'FEO', 'PG0', 'CXS',
    'AR6', 'CHD', 'WO4', 'C5P', 'UFP', 'GCP', 'HDD', 'SRT', 'STU',
    'CDP', 'TCL', '04C', 'MYA', 'URA', 'PLG', 'MTA', 'BMP', 'SAL',
    'TA1', 'UD1', 'OLA', 'BCN', 'LMR', 'BMA', 'OAA', 'TAM', 'MBO',
    'MMA', 'SPD', 'MTE', 'AP5', 'TMP', 'PGA', 'GLA', '3PG', 'FUL',
    'PQQ', '9TY', 'DUR', 'PPV', 'SPM', 'SIA', 'DUP', 'GTX', '1PG',
    'GUN', 'ETF', 'FDP', 'MFU', 'G2P', 'PC', 'DST', 'INI'
]