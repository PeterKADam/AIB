import csv
import timeit
import subprocess
import os


strings_and_scores = {
    1: ("hhppppphhppphppphp", 4),
    2: ("hphphhhppphhhhpphh", 8),
    3: ("phpphphhhphhphhhhh", 9),
    4: ("hphpphhphpphphhpphph", 9),
    5: ("hhhpphphphpphphphpph", 10),
    6: ("hhpphpphpphpphpphpphh", 9),
    7: ("pphpphhpppphhpppphhpppphh", 8),
    8: ("ppphhpphhppppphhhhhhhpphhpppphhpphpp", 14),
    9: ("pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh", 23),
    10: ("hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh", 21),
    11: ("pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp", 36),
    12: ("hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh", 42),
    13: ("hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph", 53),
    14: ("pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh", 48),
    15: ("ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh", 50),
}


PATH = "Projects/project6/"


def translate_directions(directions):
    # Define the mapping from (previous direction, current direction) to relative direction
    mapping = {
        ("n", "n"): "f",
        ("n", "e"): "r",
        ("n", "s"): "b",
        ("n", "w"): "l",
        ("e", "n"): "l",
        ("e", "e"): "f",
        ("e", "s"): "r",
        ("e", "w"): "b",
        ("s", "n"): "b",
        ("s", "e"): "l",
        ("s", "s"): "f",
        ("s", "w"): "r",
        ("w", "n"): "r",
        ("w", "e"): "b",
        ("w", "s"): "l",
        ("w", "w"): "f",
    }

    # Initialize the previous direction to north
    prev_direction = "e"

    # Translate the directions
    translated_directions = ""
    for direction in directions:
        translated_directions += mapping[(prev_direction, direction)]
        prev_direction = direction

    return translated_directions


with open(f"{PATH}/runtimes.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow(
        [
            "string number",
            "score",
            "relative score",
            "folding time",
            "absolute fold",
            "relative fold",
        ]
    )

    for key, value in strings_and_scores.items():
        data = [str(key)]
        start_time = timeit.default_timer()
        result = subprocess.run(
            f"python {PATH}/hpfold.py {value[0]}",
            shell=True,
            capture_output=True,  # Capture the output
            text=True,  # Make the output a string instead of bytes
        )
        end_time = timeit.default_timer()

        # do something here to split fold and score
        fold = result.stdout.strip()  # Get the output

        scores = subprocess.run(
            f"python {PATH}/hpview3k.py {value[0]} {''.join(fold)}", shell=True, capture_output=True, text=True
        )

        # Extract the number from the last line
        last_line = scores.stdout.strip().split("\n")[-1]

        score = int(last_line.split(": ")[-1])

        data.append(score)
        data.append(score - value[1])
        data.append(end_time - start_time)
        data.append(fold)
        data.append(translate_directions(fold))
        writer.writerow(data)
