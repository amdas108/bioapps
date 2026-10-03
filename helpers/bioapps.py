from pathlib import Path

from helpers.bio_files import (
    convert_fasta_records,
    parse_blast_queries,
    read_fasta,
    sort_descriptions,
    write_descriptions,
    write_fasta,
)

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
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Фильтрует FASTQ-последовательности по заданным параметрам.

    Аргументы:
    input_fastq: путь к входному FASTQ-файлу.
    output_fastq: имя выходного файла.
    gc_bounds: границы GC-состава в процентах.
    length_bounds: границы длины последовательности.
    quality_threshold: минимальное среднее качество.

    Возвращает:Ничего.
    """

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


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: str | None = None,
) -> None:
    """
    Переводит многострочные FASTA-последовательности в одну строку.

    Аргументы:
    input_fasta: путь к входному FASTA-файлу.
    output_fasta: путь к выходному FASTA-файлу.
     Если не указан, используется имя входного файла с суффиксом
   "_oneline".

    Возвращает:Ничего.
    """
    if output_fasta is None:
        input_path = Path(input_fasta)
        output_fasta = str(
            input_path.with_name(
                input_path.stem + "_oneline" + input_path.suffix
            )
        )

    records = read_fasta(input_fasta)
    converted_records = convert_fasta_records(records)

    write_fasta(converted_records, output_fasta)


def parse_blast_output(
    input_file: str,
    output_file: str,
) -> None:
    """
    Извлекает лучшие совпадения из результатов BLAST.

    Для каждого запроса выбирается первое название из столбца Description.
    Полученные названия сортируются по алфавиту и записываются в один столбец.

    Аргументы:
    input_file: путь к файлу с результатами BLAST.
    output_file: путь к выходному файлу.

    Возвращает:Ничего.
    s
    """
    descriptions = parse_blast_queries(input_file)
    sorted_descriptions = sort_descriptions(descriptions)

    write_descriptions(sorted_descriptions, output_file)
