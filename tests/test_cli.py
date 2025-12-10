"""Unit tests for CLI scripts and shared modules."""

import subprocess
import sys
import tempfile
from pathlib import Path

from microdna.gtf_parser import extract_basic_gene_info, extract_gene_info
from microdna.sw import sw_fill_matrix, sw_traceback


class TestGTFParser:
    """Tests for the gtf_parser module."""

    def test_extract_gene_info_complete(self):
        """All fields present should be extracted correctly."""
        attr = 'gene_id "ENSG00000123456.1"; gene_name "BRCA1"; gene_type "protein_coding"; transcript_id "ENST00000111111.1"; level 1;'
        result = extract_gene_info(attr)
        assert result["gene_id"] == "ENSG00000123456.1"
        assert result["gene_name"] == "BRCA1"
        assert result["gene_type"] == "protein_coding"
        assert result["transcript_id"] == "ENST00000111111.1"
        assert result["level"] == 1

    def test_extract_gene_info_partial(self):
        """Missing fields should return None."""
        attr = 'gene_id "ENSG00000123456.1";'
        result = extract_gene_info(attr)
        assert result["gene_id"] == "ENSG00000123456.1"
        assert result["gene_name"] is None
        assert result["gene_type"] is None

    def test_extract_gene_info_none_input(self):
        """None input should return all None values."""
        result = extract_gene_info(None)
        assert result["gene_id"] is None
        assert result["gene_name"] is None
        assert result["gene_type"] is None

    def test_extract_basic_gene_info_complete(self):
        """Basic extraction with all fields."""
        attr = 'gene_id "ENSG00000123456.1"; gene_name "BRCA1"; gene_type "protein_coding";'
        result = extract_basic_gene_info(attr)
        assert result["gene_id"] == "ENSG00000123456.1"
        assert result["gene_name"] == "BRCA1"
        assert result["gene_type"] == "protein_coding"

    def test_extract_basic_gene_info_missing(self):
        """Missing fields should return '.' (dot)."""
        attr = 'gene_id "ENSG00000123456.1";'
        result = extract_basic_gene_info(attr)
        assert result["gene_id"] == "ENSG00000123456.1"
        assert result["gene_name"] == "."
        assert result["gene_type"] == "."

    def test_extract_basic_gene_info_non_string(self):
        """Non-string input should return all dots."""
        result = extract_basic_gene_info(None)
        assert result == {"gene_id": ".", "gene_name": ".", "gene_type": "."}


class TestSmithWaterman:
    """Tests for the Smith-Waterman alignment module."""

    def test_perfect_match(self):
        """Identical sequences should have perfect alignment."""
        seq = "ATGC"
        matrix = sw_fill_matrix(seq, seq, gap=-2, mismatch=-1, match=1)
        align_a, align_b, score = sw_traceback(matrix, seq, seq, gap=-2, mismatch=-1, match=1)
        assert align_a == seq
        assert align_b == seq
        assert score == len(seq)

    def test_partial_match(self):
        """Should find local alignment in partially matching sequences."""
        seq_a = "XXATGCXX"
        seq_b = "YYATGCYY"
        matrix = sw_fill_matrix(seq_a, seq_b, gap=-2, mismatch=-1, match=1)
        align_a, align_b, score = sw_traceback(matrix, seq_a, seq_b, gap=-2, mismatch=-1, match=1)
        assert "ATGC" in align_a
        assert score >= 4  # At least the ATGC match

    def test_no_match(self):
        """Completely different sequences should have score 0 or low."""
        seq_a = "AAAA"
        seq_b = "TTTT"
        matrix = sw_fill_matrix(seq_a, seq_b, gap=-2, mismatch=-1, match=1)
        align_a, align_b, score = sw_traceback(matrix, seq_a, seq_b, gap=-2, mismatch=-1, match=1)
        assert score == 0

    def test_gap_in_alignment(self):
        """Test that gaps are correctly introduced."""
        seq_a = "ATGC"
        seq_b = "AGC"  # Missing T
        matrix = sw_fill_matrix(seq_a, seq_b, gap=-2, mismatch=-1, match=1)
        align_a, align_b, score = sw_traceback(matrix, seq_a, seq_b, gap=-2, mismatch=-1, match=1)
        # The alignment should include the common AGC subsequence
        assert score >= 2


class TestCLIScripts:
    """Tests for CLI script execution and argument parsing."""

    def test_clean_bed_help(self):
        """clean_bed.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.clean_bed", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "BED6" in result.stdout or "BED" in result.stdout

    def test_unique_eccdna_table_help(self):
        """unique_eccdna_table.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.unique_eccdna_table", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout
        assert "--output" in result.stdout

    def test_cleaner_tsv_help(self):
        """cleaner_tsv.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.cleaner_tsv", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--min-score" in result.stdout
        assert "--min-overlap-bp" in result.stdout

    def test_parse_intersected_help(self):
        """parse_intersected.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.parse_intersected", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout
        assert "--output" in result.stdout

    def test_annotate_eccdna_closest_help(self):
        """annotate_eccdna_closest.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.annotate_eccdna_closest", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout

    def test_sw_help(self):
        """sw.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.sw", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--match" in result.stdout
        assert "--gap" in result.stdout

    def test_sw_execution(self):
        """sw.py should align two sequences correctly."""
        result = subprocess.run(
            [
                sys.executable,
                "-m", "microdna.sw",
                "--A", "ATGCATGC",
                "--B", "ATGCATGC",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "Score: 8" in result.stdout

    def test_missing_required_args(self):
        """Scripts should fail gracefully when required args are missing."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.unique_eccdna_table"],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "required" in result.stderr.lower() or "error" in result.stderr.lower()

    def test_clean_bed_with_sample_data(self):
        """Test clean_bed.py with sample input data."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as infile:
            # Write sample BED-like data
            infile.write("chr1\t100\t200\t10.5\n")
            infile.write("chr1\t300\t400\t20.0\n")
            infile.write("chr2\t500\t600\t15.0\n")
            input_path = infile.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as outfile:
            output_path = outfile.name

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m", "microdna.clean_bed",
                    "-i", input_path,
                    "-o", output_path,
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            assert "Successfully formatted" in result.stdout

            # Check output file
            with open(output_path) as f:
                lines = f.readlines()
                assert len(lines) == 3
                assert "ecc_00001" in lines[0]
                assert "ecc_00002" in lines[1]
                assert "ecc_00003" in lines[2]

        finally:
            Path(input_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)

    def test_annotate_eccdna_intersected_help(self):
        """annotate_eccdna_intersected.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.annotate_eccdna_intersected", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout
        assert "--output" in result.stdout

    def test_eccdna_length_distribution_help(self):
        """eccdna_length_distribution.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.eccdna_length_distribution", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout
        assert "--output" in result.stdout

    def test_eccdna_score_distribution_help(self):
        """eccdna_score_distribution.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.eccdna_score_distribution", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout
        assert "--output" in result.stdout

    def test_get_seq_help(self):
        """get_seq.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.get_seq", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--fasta" in result.stdout
        assert "--chrom" in result.stdout

    def test_get_soft_clips_help(self):
        """get_soft_clips.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.get_soft_clips", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout

    def test_print_soft_clipped_reads_help(self):
        """print_soft_clipped_reads.py --help should work."""
        result = subprocess.run(
            [sys.executable, "-m", "microdna.print_soft_clipped_reads", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--input" in result.stdout
        assert "--num-reads" in result.stdout


class TestIntegrationWithSampleData:
    """Integration tests using synthetic sample data."""

    def test_parse_intersected_with_sample_data(self):
        """Test parse_intersected.py with synthetic bedtools intersect output."""
        # Create a synthetic 17-column bedtools intersect output
        # Columns 1-6: eccDNA BED6, Columns 7-16: Annotation BED, Column 17: overlap
        sample_data = (
            "chr1\t1000\t1200\tecc_00001\t25.5\t+\t"
            "chr1\t1050\t1500\tENSG00000001.1\t.\t+\tHAVANA\tgene\t.\t"
            'gene_id "ENSG00000001.1"; gene_name "BRCA1"; gene_type "protein_coding";\t'
            "100\n"
            "chr1\t2000\t2300\tecc_00002\t30.0\t-\t"
            "chr1\t2100\t2600\tENSG00000002.1\t.\t-\tHAVANA\tgene\t.\t"
            'gene_id "ENSG00000002.1"; gene_name "TP53"; gene_type "protein_coding";\t'
            "150\n"
        )

        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as infile:
            infile.write(sample_data)
            input_path = infile.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as outfile:
            output_path = outfile.name

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m", "microdna.parse_intersected",
                    "-i", input_path,
                    "-o", output_path,
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            assert "2 overlapping entries" in result.stdout

            # Verify output file structure
            with open(output_path) as f:
                lines = f.readlines()
                assert len(lines) == 3  # Header + 2 data rows
                assert "ecc_chr" in lines[0]  # Header
                assert "BRCA1" in lines[1] or "BRCA1" in lines[2]
                assert "TP53" in lines[1] or "TP53" in lines[2]

        finally:
            Path(input_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)

    def test_cleaner_tsv_with_sample_data(self):
        """Test cleaner_tsv.py with synthetic eccDNA annotation data."""
        # Create synthetic TSV with required columns
        header = "ecc_chr\tecc_start\tecc_end\tecc_id\tecc_length\tecc_score\tgene_id\tgene_start\tgene_end\toverlap_bp\n"
        row1 = "chr1\t1000\t1200\tecc_00001\t200\t25.0\tENSG00000001\t1050\t1500\t100\n"
        row2 = "chr1\t2000\t2300\tecc_00002\t300\t35.0\tENSG00000002\t2100\t2600\t150\n"
        row3 = "chr1\t3000\t3100\tecc_00003\t100\t5.0\tENSG00000003\t3050\t3200\t30\n"  # Low score

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as infile:
            infile.write(header)
            infile.write(row1)
            infile.write(row2)
            infile.write(row3)
            input_path = infile.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as outfile:
            output_path = outfile.name

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m", "microdna.cleaner_tsv",
                    "-i", input_path,
                    "-o", output_path,
                    "--min-score", "10.0",
                    "--min-overlap-bp", "50",
                    "--min-overlap-pct", "10.0",
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0
            assert "Cleaned data saved to" in result.stdout

            # Low score row should be filtered out
            with open(output_path) as f:
                lines = f.readlines()
                # Header + 2 rows (ecc_00001 and ecc_00002 pass filters)
                assert len(lines) == 3
                assert "ecc_00003" not in "".join(lines)

        finally:
            Path(input_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)

    def test_unique_eccdna_table_with_sample_data(self):
        """Test unique_eccdna_table.py with synthetic data."""
        # Create synthetic TSV with all required columns for aggregation
        header = "ecc_id\tecc_chr\tecc_start\tecc_end\tecc_length\tecc_score\tgene_id\tgene_name\tgene_type\n"
        row1 = "ecc_00001\tchr1\t1000\t1200\t200\t25.0\tENSG00000001\tBRCA1\tprotein_coding\n"
        row2 = "ecc_00001\tchr1\t1000\t1200\t200\t25.0\tENSG00000002\tBRCA2\tprotein_coding\n"  # Same eccDNA, different gene
        row3 = "ecc_00002\tchr1\t2000\t2300\t300\t35.0\tENSG00000003\tTP53\tprotein_coding\n"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as infile:
            infile.write(header)
            infile.write(row1)
            infile.write(row2)
            infile.write(row3)
            input_path = infile.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as outfile:
            output_path = outfile.name

        try:
            result = subprocess.run(
                [
                    sys.executable,
                    "-m", "microdna.unique_eccdna_table",
                    "-i", input_path,
                    "-o", output_path,
                ],
                capture_output=True,
                text=True,
            )
            assert result.returncode == 0

            # Should aggregate to 2 unique eccDNA
            with open(output_path) as f:
                lines = f.readlines()
                assert len(lines) == 3  # Header + 2 unique eccDNA

        finally:
            Path(input_path).unlink(missing_ok=True)
            Path(output_path).unlink(missing_ok=True)

