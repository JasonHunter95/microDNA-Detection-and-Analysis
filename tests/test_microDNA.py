"""Unit tests for microDNA.py functions."""




from microdna.microDNA import (
    most_common_char,
    count_matching_chars,
    get_consensus,
)


class TestMostCommonChar:
    """Tests for the most_common_char function."""

    def test_empty_list_returns_none(self):
        """An empty list should return None."""
        assert most_common_char([]) is None

    def test_single_char(self):
        """A single character list returns that character."""
        assert most_common_char(["A"]) == "A"

    def test_multiple_same_char(self):
        """Multiple identical characters returns that character."""
        assert most_common_char(["G", "G", "G"]) == "G"

    def test_tie_returns_first_most_common(self):
        """In a tie, Counter returns one of the most common (implementation-defined)."""
        result = most_common_char(["A", "T", "A", "T"])
        assert result in ("A", "T")

    def test_clear_winner(self):
        """Clear majority character is returned."""
        assert most_common_char(["A", "A", "A", "T", "G"]) == "A"


class TestCountMatchingChars:
    """Tests for the count_matching_chars function."""

    def test_no_matches(self):
        """No matching characters returns 0."""
        assert count_matching_chars("A", ["T", "G", "C"]) == 0

    def test_all_match(self):
        """All matching characters returns the list length."""
        assert count_matching_chars("A", ["A", "A", "A"]) == 3

    def test_partial_match(self):
        """Partial matches return correct count."""
        assert count_matching_chars("A", ["A", "T", "A", "G"]) == 2


class TestGetConsensus:
    """Tests for the get_consensus function."""

    def test_empty_list_returns_empty_string(self):
        """An empty sequence list returns an empty consensus."""
        assert get_consensus([]) == ""

    def test_single_sequence(self):
        """Single sequence doesn't meet threshold, returns empty."""
        # With default min_num_chars=15, a single sequence won't meet the threshold
        assert get_consensus(["ATGC"]) == ""

    def test_consensus_with_many_identical_sequences(self):
        """Many identical sequences should produce that sequence as consensus."""
        # 20 identical sequences should meet min_num_chars=15
        sequences = ["ATGC"] * 20
        result = get_consensus(sequences, min_num_chars=15, min_char_density=0.75)
        assert result == "ATGC"

    def test_consensus_with_variation(self):
        """Consensus should pick most common character at each position."""
        sequences = [
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "AAAA",
            "TAAA",  # Different first char
            "ATAA",  # Different second char
        ]
        result = get_consensus(sequences, min_num_chars=15, min_char_density=0.75)
        assert result == "AAAA"

    def test_custom_thresholds(self):
        """Custom thresholds should be respected."""
        sequences = ["ATG"] * 5
        # With min_num_chars=3, 5 sequences should work
        result = get_consensus(sequences, min_num_chars=3, min_char_density=0.5)
        assert result == "ATG"
