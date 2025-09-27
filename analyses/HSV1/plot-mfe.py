import math
import pickle
import subprocess
import sys

from joblib import Parallel, delayed

import ld

# Required apt packages: texlive texlive-latex-extra texlive-science
# Tex memory cap should be raised: https://tex.stackexchange.com/questions/484576/tex-capacity-exceeded-how-to-increase-main-memory-size-of-tex-live-on-a-debian

SAVETO = ld.PLOTS / 'MFE'
SAVETO.mkdir(exist_ok=True)

# Visualization Constants
ROTATION = 180 - 2  # 325 - circle size (degrees)
STEP_SIZE = 1_000
BASE_ANGLE_OFFSET = 70


def render_connection(start: int, end: int, style: str, total_bases: int, angle_adjust: float, stream):
    start_angle = 360 / (total_bases + angle_adjust) * start - BASE_ANGLE_OFFSET + ROTATION
    end_angle = 360 / (total_bases + angle_adjust) * end - BASE_ANGLE_OFFSET + ROTATION

    delta_angle = end_angle - start_angle
    if delta_angle <= 176:
        arc_radius = math.tan(delta_angle * math.pi / 360) * 10
        stream.write(
            f"\\draw[line width = 0.0mm] {style} ([shift=({end_angle}:10cm)]0,0) "
            f"arc ({end_angle + 90}:{start_angle + 270}:{arc_radius}cm); %% {start} {end} {total_bases}\n"
        )
    elif delta_angle >= 186:
        arc_radius = math.tan((360 - delta_angle) * math.pi / 360) * 10
        stream.write(
            f"\\draw[line width = 0.0mm] {style} ([shift=({start_angle}:10cm)]0,0) "
            f"arc ({start_angle + 90}:{end_angle - 90}:{arc_radius}cm); %% {start} {end} {total_bases}\n"
        )
    else:
        stream.write(f"\\draw[line width = 0.0mm] {style} ({start}.center) to ({end}.center);\n")


def job(input_file, saveto):
    with open(input_file, 'rb') as file_stream:
        _, _, _, sequence, structure, astems, zstems = pickle.load(file_stream)

    with open(saveto, 'w') as tex:
        # LaTeX Document Header
        tex.write(r"""
        \documentclass{standalone}
        \usepackage{tikz}
        \usetikzlibrary{shapes}
        \begin{document}
        
        \begin{tikzpicture}[base_style/.style={}, scale=2]
        """)

        total_bases = len(sequence)
        angle_adjust = total_bases // 9

        # Draw nodes
        for idx in range(total_bases):
            node_angle = 360 / (total_bases + angle_adjust) * idx - BASE_ANGLE_OFFSET + ROTATION
            tex.write(f"\\node [base_style] ({idx}) at ({node_angle}:10cm) {{}};\n")

            if idx % STEP_SIZE == 0 and idx != 0:
                tick = f"{idx // STEP_SIZE}kb"
                tex.write(f"\\node [scale=2] (L{idx}) at ({node_angle}:11cm) {{\\Huge {tick}}};\n")
                tex.write(f"\\node [scale=2] (T{idx}) at ({node_angle}:10.25cm) {{}};\n")
                tex.write(f"\\draw[thick] (T{idx}.center) -- ({idx}.center);\n")

            if idx > 0:
                tex.write(f"\\draw[thick] ({idx}.center) -- ({idx - 1}.center);\n")

        # Draw pair connections
        for color, pair_ranges in [('gray', astems), ('red', zstems)]:
            style = f"[{color},thick]"
            for left, right in pair_ranges:
                linds = range(left.start, left.end)
                rinds = reversed(range(right.start, right.end))

                for start, end in zip(linds, rinds):
                    assert structure[start] == '(' and structure[end] == ')'
                    render_connection(start, end, style, total_bases, angle_adjust, tex)

        # Add the total number of Z-stems to the title of the plot
        tex.write(f"\\node [scale=2,red] at (0, 10) {{\\Huge {len(zstems)} Z-stems}};\n")

        tex.write(r"""
        \end{tikzpicture}
        \end{document}
        """)

    subprocess.run(['pdflatex', saveto.name], check=True, cwd=saveto.parent, stderr=sys.stderr)


Parallel(n_jobs=-1, verbose=100)(
    delayed(job)(fold, SAVETO / (fold.name.split('.')[0] + '.tex'))
    for fold in ld.MFE.root.glob("*.scored-fold.pkl")
)
