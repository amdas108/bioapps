def read_fasta(input_fasta: str) -> list[list[str]]:
    """
    Читает FASTA-файл и разбивает его на записи.

    Аргументы:
    input_fasta: путь к входному FASTA-файлу.

    Возвращает: Список записей. Каждая запись содержит заголовок и последовательность.
    """
    records: list[list[str]] = []
    current_record: list[str] = []

    with open(input_fasta, "r") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_record:
                    records.append(current_record)

                current_record = [line]
            else:
                current_record.append(line)

        if current_record:
            records.append(current_record)

    return records


def write_fasta(records: list[list[str]], output_fasta: str) -> None:
    """
    Записывает FASTA-записи в файл.

    Аргументы:
    records: список FASTA-записей.
    output_fasta: путь к выходному FASTA-файлу.

    Возвращает:Ничего.
    """
    with open(output_fasta, "w") as file:
        for record in records:
            file.write(record[0] + "\n")
            file.write("".join(record[1:]) + "\n")


def convert_fasta_records(records: list[list[str]]) -> list[list[str]]:
    """
    Объединяет многострочные последовательности FASTA в одну строку.

    Аргументы:
    records: список FASTA-записей.

    Возвращает Список FASTA-записей с последовательностями в одну строку.
    """
    converted_records: list[list[str]] = []

    for record in records:
        header = record[0]
        sequence = "".join(record[1:])

        converted_records.append([header, sequence])

    return converted_records


def parse_blast_queries(input_file: str) -> list[str]:
    """
    Извлекает лучшие совпадения из результатов BLAST.

    Аргументы:
    input_file: путь к файлу с результатами BLAST.

    Возвращает: Список названий лучших совпадений.
    """
    descriptions: list[str] = []
    in_alignments = False
    skip_separator = False

    with open(input_file, "r") as file:
        for line in file:
            line = line.rstrip()

            if line.startswith("Sequences producing significant alignments:"):
                in_alignments = True
                skip_separator = True
                continue

            if in_alignments:
                if not line.strip():
                    in_alignments = False
                    continue

                if skip_separator:
                    skip_separator = False
                    continue

                if line.startswith("Description"):
                    continue

                description = line.strip()

                if description:
                    descriptions.append(description.split("  ")[0].strip())
                    in_alignments = False

    return descriptions


def sort_descriptions(descriptions: list[str]) -> list[str]:
    """
    Сортирует названия белков по алфавиту.

    Аргументы:
    descriptions: список названий белков.

    Возвращает:Отсортированный список названий.
    """
    return sorted(descriptions)


def write_descriptions(
    descriptions: list[str],
    output_file: str,
) -> None:
    """
    Записывает названия белков в один столбец.

    Аргументы:
    descriptions: список названий белков.
    output_file: путь к выходному файлу.

    Возвращает:Ничего.
    """
    with open(output_file, "w") as file:
        for description in descriptions:
            file.write(description + "\n")
