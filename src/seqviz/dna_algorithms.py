"""Educational DNA sequence algorithms used by the Streamlit interface.

These implementations are intentionally transparent and suitable for teaching.
They are not intended to replace validated clinical or production bioinformatics
pipelines.
"""
from __future__ import annotations

DNA_ALPHABET = set("ACGTN")


def validate_dna(sequence: str, label: str, max_length: int = 500) -> str:
    """Normalize whitespace/case and validate a DNA sequence."""
    normalized = "".join(sequence.upper().split())
    if not normalized:
        raise ValueError(f"{label} cannot be empty.")
    invalid = sorted(set(normalized) - DNA_ALPHABET)
    if invalid:
        raise ValueError(
            f"{label} contains invalid character(s): {', '.join(invalid)}. "
            "Use A, C, G, T, or N."
        )
    if len(normalized) > max_length:
        raise ValueError(
            f"{label} is too long for this interactive demo "
            f"(maximum {max_length} bases)."
        )
    return normalized


def lcs(a: str, b: str) -> tuple[str, int, list[list[int]]]:
    """Return one longest common subsequence, its length, and the DP matrix."""
    dp = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    for i, char_a in enumerate(a, 1):
        for j, char_b in enumerate(b, 1):
            if char_a == char_b:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    i, j = len(a), len(b)
    chars: list[str] = []
    while i and j:
        if a[i - 1] == b[j - 1]:
            chars.append(a[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            i -= 1
        else:
            j -= 1
    subsequence = "".join(reversed(chars))
    return subsequence, dp[-1][-1], dp


def edit_distance(
    a: str, b: str
) -> tuple[int, list[list[int]], list[tuple[str, str, str]]]:
    """Compute Levenshtein distance and recover the optimal edit operations.

    Operations are tuples of (operation, sequence-1 character, sequence-2
    character). A long dash represents a gap in that sequence.
    """
    dp = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    for i in range(len(a) + 1):
        dp[i][0] = i
    for j in range(len(b) + 1):
        dp[0][j] = j

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            substitution_cost = 0 if a[i - 1] == b[j - 1] else 1
            dp[i][j] = min(
                dp[i - 1][j] + 1,  # deletion
                dp[i][j - 1] + 1,  # insertion
                dp[i - 1][j - 1] + substitution_cost,
            )

    i, j = len(a), len(b)
    operations: list[tuple[str, str, str]] = []
    while i or j:
        if (
            i and j and a[i - 1] == b[j - 1]
            and dp[i][j] == dp[i - 1][j - 1]
        ):
            operations.append(("Match", a[i - 1], b[j - 1]))
            i -= 1
            j -= 1
        elif i and j and dp[i][j] == dp[i - 1][j - 1] + 1:
            operations.append(("Substitution", a[i - 1], b[j - 1]))
            i -= 1
            j -= 1
        elif i and dp[i][j] == dp[i - 1][j] + 1:
            operations.append(("Deletion", a[i - 1], "—"))
            i -= 1
        else:
            operations.append(("Insertion", "—", b[j - 1]))
            j -= 1

    operations.reverse()
    return dp[-1][-1], dp, operations


def smith_waterman(
    a: str, b: str, match: int = 2, mismatch: int = -1, gap: int = -2
) -> tuple[str, str, list[list[int]], int]:
    """Perform Smith–Waterman local alignment and return alignment, matrix, score."""
    matrix = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    best_score, best_pos = 0, (0, 0)

    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            diagonal = matrix[i - 1][j - 1] + (
                match if a[i - 1] == b[j - 1] else mismatch
            )
            up = matrix[i - 1][j] + gap
            left = matrix[i][j - 1] + gap
            matrix[i][j] = max(0, diagonal, up, left)
            if matrix[i][j] > best_score:
                best_score, best_pos = matrix[i][j], (i, j)

    i, j = best_pos
    aligned_a: list[str] = []
    aligned_b: list[str] = []
    while i > 0 and j > 0 and matrix[i][j] > 0:
        diagonal = matrix[i - 1][j - 1] + (
            match if a[i - 1] == b[j - 1] else mismatch
        )
        if matrix[i][j] == diagonal:
            aligned_a.append(a[i - 1])
            aligned_b.append(b[j - 1])
            i -= 1
            j -= 1
        elif matrix[i][j] == matrix[i - 1][j] + gap:
            aligned_a.append(a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(b[j - 1])
            j -= 1

    return (
        "".join(reversed(aligned_a)),
        "".join(reversed(aligned_b)),
        matrix,
        best_score,
    )


def naive_edit_distance(a: str, b: str) -> int:
    """Exponential recursive baseline, intended only for very short strings."""
    if not a:
        return len(b)
    if not b:
        return len(a)
    if a[-1] == b[-1]:
        return naive_edit_distance(a[:-1], b[:-1])
    return 1 + min(
        naive_edit_distance(a[:-1], b),
        naive_edit_distance(a, b[:-1]),
        naive_edit_distance(a[:-1], b[:-1]),
    )
