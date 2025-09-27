refexp = {
    "Z-RNA[CHM13v2]": [
        'B831009+HT-29_HSV-1_Z22-RIP_1 vs HT-29_HSV-1_Z22-input_1',
        'B831009+HT-29_HSV-1_Z22-RIP_2 vs HT-29_HSV-1_Z22-input_2',
        'B831009+HT-29_HSV-1_Z22-RIP_3 vs HT-29_HSV-1_Z22-input_3',
        'B831009+HT-29_HSV-1_Z22-RIP_4 vs HT-29_HSV-1_Z22-input_4',
        'B831009+HT-29_HSV-1_FLAG-RIP_1 vs HT-29_HSV-1_FLAG-input_1',
        'B831009+HT-29_HSV-1_FLAG-RIP_2 vs HT-29_HSV-1_FLAG-input_2',
        'B831009+HT-29_HSV-1_FLAG-RIP_3 vs HT-29_HSV-1_FLAG-input_3',
        'B831009+HT-29_HSV-1_FLAG-RIP_4 vs HT-29_HSV-1_FLAG-input_4',

        'B831009+HT-29_IAV_Z22-RIP_1 vs HT-29_IAV_Z22-input_1',
        'B831009+HT-29_IAV_Z22-RIP_2 vs HT-29_IAV_Z22-input_2',
        'B831009+HT-29_IAV_Z22-RIP_3 vs HT-29_IAV_Z22-input_3',
        'B831009+HT-29_IAV_Z22-RIP_4 vs HT-29_IAV_Z22-input_4',
        'B831009+HT-29_IAV_FLAG-RIP_1 vs HT-29_IAV_FLAG-input_1',
        'B831009+HT-29_IAV_FLAG-RIP_2 vs HT-29_IAV_FLAG-input_2',
        'B831009+HT-29_IAV_FLAG-RIP_3 vs HT-29_IAV_FLAG-input_3',
        'B831009+HT-29_IAV_FLAG-RIP_4 vs HT-29_IAV_FLAG-input_4',
    ],
    "Z-RNA[GRCm39]": [
        'B256178+MEF_IAV_Z22-RIP_1 vs MEF_IAV_input_1',
        'B256178+MEF_IAV_Z22-RIP_2 vs MEF_IAV_input_2',
        'B256178+MEF_IAV_FLAG-RIP_1 vs MEF_IAV_input_1',
        'B256178+MEF_IAV_FLAG-RIP_2 vs MEF_IAV_input_2',

        'B261790+MEF_HSV-1_Z22-RIP_1 vs MEF_HSV-1_input_1',
        'B261790+MEF_HSV-1_Z22-RIP_2 vs MEF_HSV-1_input_2',
        'B319096+MEF_HSV-1_FLAG-RIP_1 vs MEF_HSV-1_input_1',
        'B319096+MEF_HSV-1_FLAG-RIP_2 vs MEF_HSV-1_input_2',
        'B319096+MEF_HSV-1_FLAG-RIP_3 vs MEF_HSV-1_input_3',
        'B319096+MEF_HSV-1_FLAG-RIP_4 vs MEF_HSV-1_input_4',
    ],
}

inducible = {
    "Z-RNA[CHM13v2]": {
        "HSV-1": ['Input induced [HSV-1]', 'Z22 induced [HSV-1]'],
        "IAV": ['Input induced [IAV]', 'Z22 induced [IAV]'],
    },
    "Z-RNA[GRCm39]": {
        "HSV-1": ['Input induced [HSV-1]', 'Z22 induced [HSV-1]'],
        "IAV": ['Input induced [IAV]', 'Z22 induced [IAV]'],
    }
}

enriched = {
    "Z-RNA[CHM13v2]": {
        "Z22 enriched [HSV-1]": 'Z22 enriched [HSV-1]',
        "FLAG enriched [HSV-1]": 'FLAG enriched [HSV-1]',
        "Z22 enriched [IAV]": 'Z22 enriched [IAV]',
        "FLAG enriched [IAV]": 'FLAG enriched [IAV]',
    },
    "Z-RNA[GRCm39]": {
        "Z22 enriched [HSV-1]": 'Z22 enriched [HSV-1]',
        "FLAG enriched [HSV-1]": 'FLAG enriched [HSV-1]',
        "Z22 enriched [IAV]": 'Z22 enriched [IAV]',
        "FLAG enriched [IAV]": 'FLAG enriched [IAV]',
    },
}


class locations:
    palette = {
        "Elongated 3` end": "#ff8aa7",
        "De novo intergenic RNA": "#8cd7f8",
        "Intronic": "#c0c0ff",
        "Other": "#7e52a1",

        "Not significant": "#DDDDDD",
    }

    order = ["Elongated 3` end", "De novo intergenic RNA", "Intronic", "Other", "Not significant"]

    mapping = {
        'intronic': 'Intronic',
        'Intron': 'Intronic',

        'Znf-cluster[Elongated 3` end]': 'Elongated 3` end',
        'U1-cluster[Elongated 3` end]': 'Elongated 3` end',
        'Inverted histones': 'Elongated 3` end',

        "Antisense RNA": "De novo intergenic RNA",
        "Divergent transcript": "De novo intergenic RNA",
        'U1-cluster[Divergent transcript]': "De novo intergenic RNA",
        "Unclear intergenic": "De novo intergenic RNA",
        "De novo intergenic": "De novo intergenic RNA",
        'intergenic': 'De novo intergenic RNA',

        'exonic': 'Other',
        'Exon': 'Other',
        '5`UTR': 'Other',

        'Unclear': 'Other',
        'Intron-Exon': 'Other',
        '5s-cluster': 'Other',
        'tRNA-cluster': 'Other',
        'Unannotated RNA': 'Other',
        'MT': 'Other',
    }


class sequences:
    class cls:
        palette = {
            'Inv. LINE': '#f89c74',
            'Inv. SINE': '#87c55f',
            'Inv. LTR': '#dcb0f2',
            'Inv. non-repetitive': '#f6cf71',
            'Inv. other': '#66c5cc',
            'Unresolved': '#b3b3b3',
        }
        order = {
            "Z-RNA[GRCm39]": [
                'Inv. LINE', 'Inv. SINE', 'Inv. LTR', 'Inv. non-repetitive', 'Inv. other', 'Unresolved'
            ],
            "Z-RNA[CHM13v2]": [
                'Inv. SINE', 'Inv. LINE', 'Inv. non-repetitive', 'Inv. LTR', 'Inv. other', 'Unresolved'
            ],
        }


class loopsize:
    maxsize = 20_000
    binwidth = 500


class volcano:
    xymax = {
        "Z-RNA[CHM13v2]": {
            # "Z22 enriched [HSV-1]": (None, None),
            # "FLAG enriched [HSV-1]": (None, None),

            # "Z22 enriched [IAV]": (None, None),
            # "FLAG enriched [IAV]": (None, None),
            "Z22 enriched [HSV-1]": (4, 60),
            "FLAG enriched [HSV-1]": (4, 60),

            "Z22 enriched [IAV]": (4, 80),
            "FLAG enriched [IAV]": (4, 80),
        },
        "Z-RNA[GRCm39]": {
            # "Z22 enriched [HSV-1]": (None, None),
            # "FLAG enriched [HSV-1]": (None, None),

            # "Z22 enriched [IAV]": (None, None),
            # "FLAG enriched [IAV]": (None, None),
            "Z22 enriched [HSV-1]": (7, 60),
            "FLAG enriched [HSV-1]": (4, 60),

            "Z22 enriched [IAV]": (7, 60),
            "FLAG enriched [IAV]": (4, 30),
        },
    }


class thresholds:
    class induced:
        log2fc = 0.585
        padj = 0.1

    class enriched:
        log2fc = 0.585
        padj = 0.01
