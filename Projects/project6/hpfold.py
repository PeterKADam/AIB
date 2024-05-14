import sys
from typing import List, Tuple


def main(input: str):

    evenodd = [
        (
            "e"
            if i % 2 == 0 and char == "h"
            else "o" if i % 2 == 1 and char == "h" else "p"
        )
        for i, char in enumerate(input.lower())
    ]

    matches: list[Tuple[int, int]] = []
    j = len(evenodd) - 1

    matches_n = min(evenodd.count("o"), evenodd.count("e"))

    for char1 in range(
        len(evenodd),
    ):
        if evenodd[char1] == "e" and matches_n > 0:
            for char2 in range(j, 0, -1):

                if evenodd[char2] == "o":

                    matches.append((char1, char2))
                    j = char2 - 1
                    matches_n -= 1
                    break

    print(matches)
    print(len(matches))
    return


if __name__ == "__main__":
    main(sys.argv[1])
