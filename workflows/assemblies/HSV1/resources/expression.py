_expr_stages = {
    "Latency": [
        "LAT",
        "LAT-1", "LAT-2",
        "LAT-1-iso1", "LAT-2-iso1"
    ],
    "Immediate Early": [
        'RL2', 'TRL2', 'IRL2',
        'RS1', 'IRS1', 'TRS1',
        'UL54', 'US1', 'US12'
    ],
    "Early": [
        'UL2', 'UL5', 'UL8', 'UL9', 'UL12', 'UL23', 'UL29', 'UL30', 'UL39',
        'UL40', 'UL42', 'UL50', 'UL52', 'US3'
    ],
    "Late": [
        # https://doi.org/10.1128/JVI.76.16.8003-8010.2002
        "AL-RNA",
        # https://www.microbiologyresearch.org/content/journal/jgv/10.1099/0022-1317-79-12-3033
        'UL4',
        # https://pubmed.ncbi.nlm.nih.gov/7831802/
        'UL6', 'UL7',
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2224439/
        'UL14',
        # https://pubmed.ncbi.nlm.nih.gov/10639317/
        'UL33',
        # https://www.sciencedirect.com/science/article/pii/S1535947620337385
        'UL43',
        # https://www.mcponline.org/article/S1535-9476(20)33738-5/fulltext#supplementaryMaterial
        'UL55', 'UL56',

        'TRL1', 'IRL1', 'RL1',
        'UL1', 'UL3', 'UL10', 'UL11', 'UL13', 'UL15', 'UL16', 'UL17', 'UL18', 'UL19', 'UL20', 'UL21', 'UL22', 'UL24',
        'UL25', 'UL26', 'UL26.5', 'UL27', 'UL28', 'UL31', 'UL32', 'UL34', 'UL35', 'UL36', 'UL37', 'UL38', 'UL41',
        'UL44', 'UL45', 'UL46', 'UL47', 'UL48', 'UL49', 'UL49A', 'UL51', 'UL53',

        'US2', 'US4', 'US5', 'US6', 'US7', 'US8', 'US8A', 'US9', 'US10', 'US11',

        # Clusters
        "UL19/UL20", "UL33/UL34", "US10/US11", "UL13/UL14"
    ]
}
rev = {}
for group, genes in _expr_stages.items():
    for g in genes:
        assert g not in rev, g
        rev[g] = group
STAGES = rev
