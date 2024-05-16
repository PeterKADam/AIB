import sys
from typing import List, Tuple


def hp_matches(input: str):

    evenodd = [
        ("e" if i % 2 == 0 and char == "h" else "o" if i % 2 == 1 and char == "h" else "p")
        for i, char in enumerate(input.lower())
    ]

    matches: list[Tuple[int, int]] = []
    j = len(evenodd) - 1

    matches_n = min(evenodd.count("o"), evenodd.count("e"))

    for char1 in range(
        len(evenodd),
    ):
        if evenodd[char1] == "o" and matches_n > 0:
            for char2 in range(j, 0, -1):

                if evenodd[char2] == "e" and char1 < char2 and (char2 - char1) > 1:

                    matches.append((char1, char2))
                    j = char2 - 1
                    matches_n -= 1
                    break

    return matches


def fold(string):
    length = len(string)
    match = hp_matches(string)
    folding = []

    if len(match) == 1:
        folding.append(match[0][0] * "e")
        center = match[0][1] - match[0][0]
        east = (center - 1) // 2
        folding.append(east * "e")
        folding.append("s")
        folding.append(east * "w")

    else:
        folding.append((match[0][0]) * "e")

        for i in range(1, len(match)):
            ps = match[i][0] - match[i - 1][0]
            if ps >= 3:
                north = (ps // 2) - 1
                folding.append(f"{north * 'n'}e{north * 's'}e")
            else:
                folding.append(ps * "e")
            if i == len(match) - 1:
                center = match[i][1] - match[i][0]
                if center >= 3:
                    east = (center - 1) // 2
                    folding.append(f"{east * 'e'}s{east * 'w'}")
                else:
                    folding.append("esw")

        for j in range(len(match) - 1, 0, -1):
            ps = match[j - 1][1] - match[j][1]
            if ps >= 3:
                south = (ps - 1) // 2
                folding.append(f"{south * 's'}w{south * 'n'}w")
            else:
                folding.append(ps * "w")

    end = length - match[0][1] - 1
    folding.append(end * "w")

    # print(len("".join(folding)))
    return "".join(folding)


def main(input: str):
    print(fold(input))


if __name__ == "__main__":
    main(sys.argv[1])
