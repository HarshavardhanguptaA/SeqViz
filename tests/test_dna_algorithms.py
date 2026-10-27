"""Unit tests for the educational DNA alignment algorithms."""
import pytest

from seqviz.dna_algorithms import (
    edit_distance,
    lcs,
    naive_edit_distance,
    smith_waterman,
    validate_dna,
)


def test_validate_dna_normalizes_case_and_whitespace():
    assert validate_dna(" ac gt\nN ", "Sequence") == "ACGTN"


@pytest.mark.parametrize("sequence", ["", "   ", "ACGX", "123"])
def test_validate_dna_rejects_empty_or_invalid_input(sequence):
    with pytest.raises(ValueError):
        validate_dna(sequence, "Sequence")


def test_validate_dna_enforces_length_limit():
    with pytest.raises(ValueError, match="maximum 3"):
        validate_dna("ACGT", "Sequence", max_length=3)


def test_lcs_finds_common_subsequence():
    subsequence, length, matrix = lcs("ABCBDAB", "BDCABA")
    assert length == 4
    assert len(subsequence) == length
    assert _is_subsequence(subsequence, "ABCBDAB")
    assert _is_subsequence(subsequence, "BDCABA")
    assert len(matrix) == 8
    assert len(matrix[0]) == 7


def test_edit_distance_exact_match():
    distance, matrix, operations = edit_distance("ACGT", "ACGT")
    assert distance == 0
    assert all(operation[0] == "Match" for operation in operations)
    assert len(matrix) == 5 and len(matrix[0]) == 5


@pytest.mark.parametrize(
    ("a", "b", "expected"),
    [
        ("kitten", "sitting", 3),
        ("GATTACA", "GATACA", 1),
        ("ACGT", "ACGGT", 1),
        ("A", "T", 1),
        ("", "AC", 2),
        ("AC", "", 2),
    ],
)
def test_edit_distance_mutation_cases(a, b, expected):
    assert edit_distance(a, b)[0] == expected


def test_edit_operations_include_insert_delete_and_substitution():
    assert edit_distance("AC", "ATC")[0] == 1
    assert edit_distance("ATC", "AC")[0] == 1
    assert edit_distance("AC", "AT")[0] == 1


def test_smith_waterman_returns_consistent_alignment():
    aligned_a, aligned_b, matrix, score = smith_waterman(
        "TTACGTAA", "GGACGTCC"
    )
    assert aligned_a == aligned_b == "ACGT"
    assert score == 8
    assert len(aligned_a) == len(aligned_b)
    assert len(matrix) == 9 and len(matrix[0]) == 9


def test_naive_and_dynamic_programming_agree_for_short_inputs():
    examples = [("", ""), ("A", "T"), ("AC", "AG"), ("GAT", "GCT")]
    for a, b in examples:
        assert naive_edit_distance(a, b) == edit_distance(a, b)[0]


def _is_subsequence(candidate, sequence):
    iterator = iter(sequence)
    return all(any(char == current for current in iterator) for char in candidate)
