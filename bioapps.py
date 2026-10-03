from pathlib import Path

from helpers.fastq import (
    check_gc,
    check_length,
    check_quality,
    read_fastq_record,
    validate_bounds,
    validate_quality_threshold,
    write_fastq_record,
)


def filter_fastq(
    input_fastq,
    output_fastq,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
    validate_bounds(gc_bounds)
    validate_bounds(length_bounds)
    validate_quality_threshold(quality_threshold)

    output_dir = Path("filtered")
    output_dir.mkdir(exist_ok=True)

    output_path = output_dir / output_fastq

    with open(input_fastq, "r") as input_file, \
            open(output_path, "w") as output_file:

        while True:
            record = read_fastq_record(input_file)

            if record is None:
                break

            header, sequence, plus, quality = record

            if not check_gc(sequence, gc_bounds):
                continue

            if not check_length(sequence, length_bounds):
                continue

            if not check_quality(quality, quality_threshold):
                continue

            write_fastq_record(output_file, record)


