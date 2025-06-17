import tempfile
from pathlib import Path
from click.testing import CliRunner

from hicue.cli.extract import extract
from hicue.cli.extract2d import extract2d
from hicue.cli.tracks import tracks
from hicue.cli.compare import compare


def test_hicue_extract():
    """Test the hicue extract CLI command with actual test data."""
    runner = CliRunner()

    # Get the project root directory
    project_root = Path(__file__).parent.parent
    tests_data_dir = project_root / "tests_data"
    bed_file = tests_data_dir / "microC_Constantino_anchors_subset.bed"
    cool_file = tests_data_dir / "microC_Constantino.cool"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "out"

        # Run the hicue extract command
        result = runner.invoke(
            extract,
            [
                str(output_dir),
                str(bed_file),
                str(cool_file),
                "--pileup",
            ],
        )

        # Check that the command executed successfully
        assert (
            result.exit_code == 0
        ), f"Command failed with exit code {result.exit_code}. Output: {result.output}"

        # Check that output directory was created
        assert output_dir.exists(), "Output directory should have been created"

        # Check that some output files were generated
        output_files = list(output_dir.rglob("*"))
        assert len(output_files) > 0, "No output files were generated"

        # Check for expected output structure (log file should exist)
        log_files = list(output_dir.glob("*_log.txt"))
        assert len(log_files) > 0, "Log file should have been created"

        print(
            f"Test passed! Generated {len(output_files)} output files in {output_dir}"
        )


def test_hicue_extract_loops():
    """Test the hicue extract CLI command with actual test data."""
    runner = CliRunner()

    # Get the project root directory
    project_root = Path(__file__).parent.parent
    tests_data_dir = project_root / "tests_data"
    bed_file = tests_data_dir / "microC_Constantino_anchors_subset.bed"
    cool_file = tests_data_dir / "microC_Constantino.cool"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "out"

        # Run the hicue extract command
        result = runner.invoke(
            extract,
            [
                str(output_dir),
                str(bed_file),
                str(cool_file),
                "--pileup",
                "--pileup_method",
                "mean",
                "--nb_pos",
                "3",
                "--detrending",
                "patch",
                "--diag_mask",
                "5000",
                "--cmap_limits",
                "-1",
                "1",
            ],
        )

        # Check that the command executed successfully
        assert (
            result.exit_code == 0
        ), f"Command failed with exit code {result.exit_code}. Output: {result.output}"

        # Check that output directory was created
        assert output_dir.exists(), "Output directory should have been created"

        # Check that some output files were generated
        output_files = list(output_dir.rglob("*"))
        assert len(output_files) > 0, "No output files were generated"

        # Check for expected output structure (log file should exist)
        log_files = list(output_dir.glob("*_log.txt"))
        assert len(log_files) > 0, "Log file should have been created"

        print(
            f"Test passed! Generated {len(output_files)} output files in {output_dir}"
        )


def test_hicue_extract2():
    """Test the hicue extract CLI command with actual test data."""
    runner = CliRunner()

    # Get the project root directory
    project_root = Path(__file__).parent.parent
    tests_data_dir = project_root / "tests_data"
    bedpe_file = tests_data_dir / "microC_Constantino_loops_subset.bedpe"
    cool_file = tests_data_dir / "microC_Constantino.cool"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "out"

        # Run the hicue extract command
        result = runner.invoke(
            extract2d,
            [
                str(output_dir),
                str(bedpe_file),
                str(cool_file),
                "--pileup",
                "--detrending",
                "ps",
            ],
        )

        # Check that the command executed successfully
        assert (
            result.exit_code == 0
        ), f"Command failed with exit code {result.exit_code}. Output: {result.output}"

        # Check that output directory was created
        assert output_dir.exists(), "Output directory should have been created"

        # Check that some output files were generated
        output_files = list(output_dir.rglob("*"))
        assert len(output_files) > 0, "No output files were generated"

        # Check for expected output structure (log file should exist)
        log_files = list(output_dir.glob("*_log.txt"))
        assert len(log_files) > 0, "Log file should have been created"

        print(
            f"Test passed! Generated {len(output_files)} output files in {output_dir}"
        )


def test_hicue_track():
    """Test the hicue tracks CLI command with a track."""
    runner = CliRunner()

    # Get the project root directory
    project_root = Path(__file__).parent.parent
    tests_data_dir = project_root / "tests_data"
    track_file = tests_data_dir / "Scc1.bw"
    cool_file = tests_data_dir / "microC_Constantino.cool"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "out"

        # Run the hicue extract command
        result = runner.invoke(
            tracks,
            [
                "--threshold",
                "min",
                "5",
                "--pileup",
                "--pileup_method",
                "median",
                "--detrending",
                "patch",
                "--loops",
                str(output_dir),
                str(track_file),
                str(cool_file),
            ],
        )

        # Check that the command executed successfully
        assert (
            result.exit_code == 0
        ), f"Command failed with exit code {result.exit_code}. Output: {result.output}"

        # Check that output directory was created
        assert output_dir.exists(), "Output directory should have been created"

        # Check that some output files were generated
        output_files = list(output_dir.rglob("*"))
        assert len(output_files) > 0, "No output files were generated"

        # Check for expected output structure (log file should exist)
        log_files = list(output_dir.glob("*_log.txt"))
        assert len(log_files) > 0, "Log file should have been created"

        print(
            f"Test passed! Generated {len(output_files)} output files in {output_dir}"
        )


def test_hicue_compare():
    """Test the hicue compare CLI command with actual test data."""
    runner = CliRunner()

    # Get the project root directory
    project_root = Path(__file__).parent.parent
    tests_data_dir = project_root / "tests_data"
    bedpe_file = tests_data_dir / "microC_Constantino_loops_subset.bedpe"
    cool_file_1 = tests_data_dir / "microC_Constantino.cool"
    cool_file_2 = tests_data_dir / "microC_Constantino.cool"

    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir) / "out"

        # Run the hicue extract command
        result = runner.invoke(
            compare,
            [
                "--pileup",
                "--detrending",
                "ps",
                str(output_dir),
                str(bedpe_file),
                f"{str(cool_file_1)},{str(cool_file_2)}",
            ],
        )

        # Check that the command executed successfully
        assert (
            result.exit_code == 0
        ), f"Command failed with exit code {result.exit_code}. Output: {result.output}"

        # Check that output directory was created
        assert output_dir.exists(), "Output directory should have been created"

        # Check that some output files were generated
        output_files = list(output_dir.rglob("*"))
        assert len(output_files) > 0, "No output files were generated"

        # Check for expected output structure (log file should exist)
        log_files = list(output_dir.glob("*_log.txt"))
        assert len(log_files) > 0, "Log file should have been created"

        print(
            f"Test passed! Generated {len(output_files)} output files in {output_dir}"
        )
